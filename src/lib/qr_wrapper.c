/**
 * Copyright (C) Till Blaha 2022-2023
 * MAVLab -- Faculty of Aerospace Engineering -- Delft University of Techology
 */

/**
 * @file qr_wrapper.c
 * 
 * @brief Implementations of qr_wrapper.h
*/

#include "qr_wrapper.h"
#include "qr_solve/qr_solve.h"
#include <math.h>
#include <string.h>

#include <stdio.h>

//#include "f2c.h"
//#include "clapack.h"

void qr_wrapper(int m, int n, int* perm, num_t** A, num_t** Q, num_t** R) {
  num_t in[AS_N_C*AS_N_C]; // changed from VLA...
  int k = 0;
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < m; i++) {
      in[k++] = A[i][perm[j]];
    }
  }
  int jpvt[AS_N_U];
  //int kr;
  num_t tau[AS_N_U];
  num_t work[AS_N_U];
  int job = 0;

  dqrdc ( in, m, m, n, tau, jpvt, work, job );
  for (int i = 0; i < m*n; i++) {
    if (i%m > i/m) {
      R[i%m][i/m] = 0;
    } else {
      R[i%m][i/m] = in[i];
    }
  }
  //sgeqr2_ ( &m, &n, in, &m, tau, work, &job );

  //num_t Qout[AS_N_C*AS_N_C];
  //dorgqr ( m, n, in, Qout, tau);
  //org2r ( m, m, n, in, tau, work, &job );

  #ifdef debug_qr
  printf("in = [\n");
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      printf("%f ", in[i + j*m]);
    }
    printf(";\n");
  }
  printf("];\n\n");

  printf("tau = [");
  for (int i = 0; i < m; i++) {
    printf("%f ", tau[i]);
    printf(";\n");
  }
  printf("];\n\n");
  #endif

  org2r ( m, n, in, tau );

  #ifdef debug_qr
  printf("in = [\n");
  for (int i = 0; i < m; i++) {
    printf("[");
    for (int j = 0; j < m; j++) {
      printf("%f ,", in[i + j*m]);
    }
    printf("],\n");
  }
  printf("];\n\n");

  #endif

  for (int i = 0; i < m*m; i++)
    //Q[i%m][i/m] = Qout[i];
    Q[i%m][i/m] = in[i];

}

void qr_wrapper2(int m, int n, num_t* in, num_t* R) {
  int jpvt[AS_N_U];
  //int kr;
  num_t tau[AS_N_U];
  num_t work[AS_N_U];
  int job = 0;

  dqrdc ( in, m, m, n, tau, jpvt, work, job );

  // Copy the upper triangular part of in to R
  for (int j = 0; j < n; j++)  // Loop over columns
      for (int i = 0; i <= j && i < m; i++)  // Loop over rows, only up to the diagonal (i <= j)
          R[j * m + i] = in[j * m + i];

  #ifdef debug_qr
  printf("in = [\n");
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      printf("%f ", in[i + j*m]);
    }
    printf(";\n");
  }
  printf("];\n\n");

  printf("tau = [");
  for (int i = 0; i < m; i++) {
    printf("%f ", tau[i]);
    printf(";\n");
  }
  printf("];\n\n");
  #endif

  org2r ( m, n, in, tau );

  #ifdef debug_qr
  printf("in = [\n");
  for (int i = 0; i < m; i++) {
    printf("[");
    for (int j = 0; j < m; j++) {
      printf("%f ,", in[i + j*m]);
    }
    printf("],\n");
  }
  printf("];\n\n");

  #endif

}

int dorgqr ( int m, int n, const num_t a[], num_t q[], num_t tau[])
{
  // initialize to 0
  // TODO: is this necessary?
  for (int i = 0; i < m*m; i++) {
    q[i] = 0.;
  }
  for (int k=0; k<m; k++) {
    // start with identity
    q[k + k*m] = 1.0F;
  }
  num_t sqtau, isqtau, tsum;
  for (int i=0; i<n; i++) {
    //sqtau = (tau[i] > 0.f) ? sqrtf(tau[i]) : 0.;
    //isqtau = (sqtau > 1e-8f) ? 1.f/sqtau : 0.;
    sqtau = sqrtf(tau[i]); //tau guaranteed to be non negative
    isqtau = 1/sqtau;

    for (int k=0; k<m; k++) {
      // tsum = Q(k, :)*v
      tsum = q[k + i*m] * sqtau; // *v(i), but v(i)=1
      for (int l=i+1; l<m; l++) {
        tsum += q[k + l*m] * a[l + i*m] * isqtau;
      }

      // Q(k, j) = Q(k, j) - tau(i)*(Q(k, :)*v)*v(j)
      q[k + i*m] -= tsum * sqtau; // *v(j), but v(j==i)=1
      for (int j=i+1; j<m; j++) {
        q[k + j*m] -= tsum*a[j + i*m] * isqtau;
      }
    }
  }
  // i'll eat my hat if this works first try
  return 0;
}

int org2r ( int m, int k, num_t* A, num_t* TAU )
{
/*  Purpose */
/*  ======= */

/*  SORG2R generates an m by n real matrix Q with orthonormal columns, */
/*  which is defined as the first n columns of a product of k elementary */
/*  reflectors of order m */

/*        Q  =  H(1) H(2) . . . H(k) */

/*  as returned by SGEQRF. */

/* Till Blaha 2024 */

    // NOTE: full Q will be generated, so must be M x M

    if (k > m) return 1;

    int ik, in, im;
    int lda = m;
    num_t tmp;

    for ( ik = 0; ik < k; ik++ ) A[ik + ik*lda] = TAU[ik];
    for ( im = k; im < m; im++ ) A[im + im*lda] = 1.0f;
    for ( ik = 0; ik < k; ik++ ) TAU[ik] = 1. / (TAU[ik]);

    // last householder factor k-1
    ik = k-1;
    // Q <-- I - TAU[k-1] * outer(vk-1, vk-1) * I
    for ( im = m-1; im >= ik; im-- ) {
        tmp = A[im + ik*lda];
        for ( in = ik; in < m; in++ )
            A[in + im*lda] = ((float) (in == im)) - TAU[ik] * A[in + ik*lda] * tmp;
    }

    ik--;
    for ( ; ik >= 0; ik-- ) {
        // Q <-- Q - TAU[ik] * outer(v[k], v[k]) * Q
        // note, that only lower triangle of Q[ik-1:, ik-1:] is actually Q at the
        // start of each iteration. The rest still contains the householder
        // vectors 0:ik-1 on the lower triangular (including diagonal) and the 
        // upper triangular (excluding diagonal) still contains R
        for ( im = m-1; im >= ik; im-- ) {
            // inner product v[ik]**T Q[:, im] will be stored in tmp
            if ( ik == im )
                tmp = A[ik + ik*lda];
            else {
                tmp = A[ik+1 + ik*lda] * A[ik+1 + im*lda];
                for ( in = ik+2; in < m; in++ )
                    tmp += A[in + ik*lda] * A[in + im*lda];
            }
            tmp *= TAU[ik];
            for ( in = ik; in < m; in++ ) {
                if ( ik == im )
                    A[in + im*lda] = ((float)(in == im)) - A[in + ik*lda] * tmp;
                else {
                    if ( in == ik )
                        A[in + im*lda] = -A[in + ik*lda] * tmp;
                    else
                        A[in + im*lda] = A[in + im*lda] - A[in + ik*lda] * tmp;
                }
            }
        }
    }

    return 0;
}


