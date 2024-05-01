/**
 * Copyright (C) Till Blaha 2022-2023
 * MAVLab -- Faculty of Aerospace Engineering -- Delft University of Techology
 */

/**
 * @file setupWLS.c
 * 
 * @brief Implementations of setupWLS.h
 */

#include "setupWLS.h"
#include <math.h>

void setupWLS_A(
    const num_t B[AS_N_V*AS_N_U], const num_t Wv[AS_N_V], num_t Wu[AS_N_U],
    const int n_v, const int n_u,
    const num_t theta, const num_t cond_bound,
    num_t A[AS_N_C*AS_N_U], num_t* gamma){

    int n_c = n_v + n_u;

    // calculate (B * Wv * Wv * B^T) exploting the symmetry
    num_t A2[AS_N_V][AS_N_V];
    num_t * A2_ptr[AS_N_V];
    for (int i=0; i<n_v; i++)
        A2_ptr[i] = A2[i];

    for (int i=0; i<n_v; i++) {
        for (int j=i; j<n_v; j++) {
            A2_ptr[i][j] = 0;
            for (int k=0; k<n_u; k++)
                A2_ptr[i][j] += B[i+k*n_v]*B[j+k*n_v];

            A2_ptr[i][j] *= Wv[i]*Wv[i];
            if (i != j) {
              A2_ptr[j][i] = A2_ptr[i][j];
            }
        }
    }

    // get minimum of Wu and normalize by that. Also get max of resultant
    num_t min_diag = +INFINITY;
    num_t max_diag = 0.;
    for (int i=0; i<n_u; i++) {
        min_diag = (min_diag > Wu[i]) ? Wu[i] : min_diag;
        max_diag = (max_diag < Wu[i]) ? Wu[i] : max_diag;
    }

    min_diag = (min_diag > 1e-6) ? min_diag : 1e-6;
    num_t min_diag_inv = 1/min_diag;
    for (int i=0; i<n_u; i++)
        Wu[i] *= min_diag_inv;

    max_diag *= min_diag_inv; 
    // has to be 1 or larger than 1, so divisions by max_diag are fine

    num_t max_sig;
    if (cond_bound > 0) {
        // estimate minimum gamma, and then use it, unless theta results in a 
        // higher value
        gamma_estimator(n_v, A2_ptr, cond_bound, gamma, &max_sig);
        *gamma = (*gamma > sqrtf(max_sig)*theta/max_diag) ? *gamma 
        : sqrtf(max_sig)*theta/max_diag;
    } else {
        // bypass the gamma returned from the estimator, but use the max_sig
        num_t dummy;
        gamma_estimator(n_v, A2_ptr, 1., &dummy, &max_sig);
        *gamma = sqrtf(max_sig)*theta/max_diag;
    }

    // add sparse part to A (from actuator penalty)
    for (int i=0; i<n_c; i++) {
        if (i < n_v) {
            for (int j=0; j<n_u; j++)
                A[i+j*n_c] = Wv[i]*B[i+j*n_v];
        } else {
            for (int j=0; j<n_u; j++) {
                if ( (i-n_v) == j ) {
                    A[i+j*n_c] = *gamma*Wu[i-n_v];
                } else {
                    A[i+j*n_c] = 0;
                }
            }
        }
    }
}

void setupWLS_b(
    const num_t v[AS_N_V], const num_t ud[AS_N_U],
    const num_t Wv[AS_N_V], const num_t Wu_norm[AS_N_U],
    const int n_v, const int n_u,
    const num_t gamma,
    num_t b[AS_N_C]){

    int n_c = n_v + n_u;

    for (int i=0; i<n_c; i++) {
        if (i < n_v) {
            b[i] = Wv[i]*v[i];
        } else {
            b[i] = gamma*Wu_norm[i-n_v]*ud[i-n_v];
        }
    }
}

void genAbFromCtlAlloc(
    int n_v, int n_u,
    num_t G[AS_N_V*AS_N_U], num_t Wv[AS_N_V], num_t Wu[AS_N_U], num_t up[AS_N_U],
    num_t dv[AS_N_V],
    num_t A[AS_N_C*AS_N_U], num_t b[AS_N_C])
{
    num_t theta = 2.0e-9;
    num_t cond_bound = 4e7;
    num_t gamma;
    setupWLS_A(G, Wv, Wu, n_v, n_u, theta, cond_bound, A, &gamma);
    setupWLS_b(dv, up, Wv, Wu, n_v, n_u, gamma, b);
}

void genHbetaFromAb(
    int m, int n, num_t A[AS_N_C*AS_N_U], num_t b[AS_N_C],
    num_t H[AS_N_U*AS_N_U], num_t beta[AS_N_U])
{
    // H = 0.5 A**T * A
    // beta = -A**T * b

    int row, col, k;
    num_t tmp;
    for (row = 0; row < n; row++) {
        for (col = 0; col < n; col++) {
            tmp = 0.;
            for (k = 0; k < m; k++)
                tmp += A[row*m + k] * A[col*m + k];
            H[col*n + row] = tmp;
        }
    }

    for (row = 0; row < n; row++) {
        tmp = 0.;
        for (k = 0; k < m; k++)
            tmp -= A[row*m + k] * b[k];
        beta[row] = tmp;
    }
}

void gamma_estimator(
    const int n, num_t** A2, const num_t cond_target,
    num_t* gamma, num_t* max_sig){
    /*
    % returns gamma to meet upper bound on condition number. Should run in
    % O(n*d^2) time when optimised. Also return upper bound on maximum
    % eigenvalue
    */
    *max_sig = 0;

    num_t R;
    for (int i=0; i<n; i++) {
        R = 0;
        for (int j=0; j<n; j++) {
            if (j != i)
#ifdef AS_SINGLE_FLOAT
                R += fabsf(A2[i][j]);
#else
                R += fabs(A2[i][j]);
#endif
        }
        if (*max_sig < (A2[i][i]+R))
            *max_sig = A2[i][i]+R;
    }

    *gamma = sqrtf(*max_sig / cond_target);
}

void cond_estimator(
    const int d, num_t** A2, const num_t min_diag2,
    num_t* cond_est, num_t* max_sig){
    /*
    % returns upper bound on condition number of A matrix. Should run in
    % O(n*d^2) time when optimised
    */
    *max_sig = 0;

    num_t R;
    for (int i=0; i<d; i++) {
        R = 0;
        for (int j=0; j<d; j++) {
            if (j != i)
#ifdef AS_SINGLE_FLOAT
                R += fabsf(A2[i][j]);
#else
                R += fabs(A2[i][j]);
#endif
        }
        *max_sig = (*max_sig < A2[i][i]+R) ? A2[i][i]+R : *max_sig;
    }

    *cond_est = *max_sig / ((min_diag2 > 1e-10) ? min_diag2 : 1e-10);

}