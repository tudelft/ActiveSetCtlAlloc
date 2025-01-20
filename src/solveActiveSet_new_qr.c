/** @file solveActiveSet_qr.c
 * @brief Use active-set QR algorithm to solve constrained least squares
 *
 * This algorithm will find the optimal solution to the least squares problem 
 * setup with setup_wls.c
 *
 * Uses an active-set algorithm with QR subspace optimisation. The QR components
 * are updated at every step instead of naively re-factorised.
 *
 * written by Till Blaha
 * MAVLab Delft University of Technology
 */

#include "solveActiveSet.h"
#include <stdio.h>
/*#include "std.h"*/
#include <inttypes.h>
#include <stdbool.h>
#include <math.h>

#include <string.h>
#include <math.h>
#include <float.h>
#include "qr_wrapper.h"
#include "qr_updates.h"
#include "qr_solve/qr_solve.h"
#include "qr_solve/r8lib_min.h"
#include "sparse_math.h"


/**
 * @brief active set algorithm for constrained least squares, QR factorisation
 *
 * Takes a representation of the problem
 * 
 *    min || Au - b ||^2 
 *      st. lb <= u <= lb
 * 
 * and returns the optimal solution u and the active-set, while respective a
 * certain iteration limit.
 * 
 * Inputs to this function can be derived from setup_wls.c
 * 
 * @param A column major representation of A
 * @param b the rhs
 * @param lb The minimum u vector
 * @param ub The maximum u vector
 * @param x Input/Output: initial guess to be iteratively improved
 * @param Ws Input/Output: indicates whether x[i] is bounded at lb[i] (-1),
 * at ub[i] (+1), or not bounded (0)
 * @param updating unused. Could be used to switch between updating and naive 
 * factorisations
 * @param imax Maximum number of iterations taken by the algorithm 
 * @param n Length of u
 * @param n_v Lenght of v
 * @param iter Iterations taken
 * 
 * @note n+n_v must be <= 256
 *
 * @return 0: exact solution found. 1: imax hit, solution not optimal
 */
activeSetExitCode solveActiveSet_new_qr(
  const num_t A_col[AS_N_C*AS_N_U], const num_t b[AS_N_C],
  const num_t lb[AS_N_U], const num_t ub[AS_N_U], num_t x[AS_N_U],
  int8_t Ws[AS_N_U], int imax, const int n_u, const int n_v,
  int *iter, int *n_free, num_t costs[])
{

#ifndef AS_RECORD_COST
  (void) costs;
#endif

  if(!imax) imax = 100;

  int exit_code = AS_ITER_LIMIT;

  int i, j;
  int n = n_u;
  int m = n_u + n_v;
  //int n_v = m-n;

  for (i = 0; i < n; i++)
    x[i] = (Ws[i] > 0) ? ub[i] : (Ws[i] < 0) ? lb[i] : 0.5*(lb[i] + ub[i]);

  int permutation[AS_N_U]; memset(permutation, 0, sizeof(int)*n);
  (*n_free) = 0;
  int i_bnd = 0;
  for (i = 0; i < n; i++)
    if (Ws[i] == 0) permutation[(*n_free)++] = i;

  for (i = 0; i < n; i++)
    if (Ws[i] != 0) permutation[(*n_free)+(i_bnd++)] = i;

  num_t Q[AS_N_C*AS_N_C];
  num_t R[AS_N_C*AS_N_U];
  for (i=0; i < n; i++)
    memcpy(Q+i*m, A_col+permutation[i]*m, sizeof(num_t)*m);

  // initial factorisation
  qr_wrapper2(m, n, Q, R);
  #ifdef debug_qr
  print_debug(A_ptr, Q_ptr, R_ptr, &n, &m);
  #endif

  num_t q[AS_N_U];
  num_t z[AS_N_U];
  bool nan_found = false;

  // -------------- Start loop ------------
  *iter = 0;
  while (++(*iter) <= imax) {
    num_t c[AS_N_U];
    for (i=0; i < (*n_free); i++) {
      c[i] = 0;
      for (j=0; j < m; j++) {
        c[i] += Q[j + m*i]*b[j];
      }
    }

    num_t u_bound_perm[AS_N_U];
    for (i=0; i < n - (*n_free); i++) {
      if (Ws[permutation[i+(*n_free)]] > 0)
        u_bound_perm[i] = ub[permutation[i+(*n_free)]];
      else if (Ws[permutation[i+(*n_free)]] < 0) 
        u_bound_perm[i] = lb[permutation[i+(*n_free)]];
#ifdef debug_qr
      else
        printf("W out of bounds for bounded variables");
#endif
    }

    for (i=0; i < (*n_free); i++) {
      for (j=0; j < n - (*n_free); j++) {
        c[i] -= R[i + m*((*n_free)+j)] * u_bound_perm[j];
      }
    }

    backward_tri_solve2((*n_free), m, R, c, q);
    
    for (i = 0; i < (*n_free); i++) {
      if (q[i] != q[i]) {
      //if (isnan(q[i])) {
        // break immediately with error
        nan_found = true;
        break;
      }
      z[permutation[i]] = q[i];
    }
    if (nan_found) {
      exit_code = AS_NAN_FOUND_Q;
      break;
    }
    for (i = (*n_free); i < n; i++) {
      z[permutation[i]] = x[permutation[i]];
    }

    uint8_t n_violated = 0;
    int8_t W_temp[AS_N_U];
    n_violated = check_limits_tol((*n_free), AS_CONSTR_TOL, z, lb, ub, W_temp, permutation);

    if (!n_violated) {
      for (i = 0; i < (*n_free); i++)
        x[permutation[i]] = z[permutation[i]];

      if ((*n_free) == n) {
        // no active constraints, we are optinal and feasible
#ifdef AS_RECORD_COST
        if ((*iter) <= AS_RECORD_COST_N)
          costs[(*iter)-1] = calc_cost(A, b, us, n, n_v);
#endif
        exit_code = AS_SUCCESS;
        break;
      } else {
        // active constraints, check for optimality
        num_t d[AS_N_U];
        for (i=(*n_free); i < n; i++) {
          d[i] = 0;
          for (j=0; j < m; j++) {
            d[i] += Q[j + m*i]*b[j];
          }
        }

        for (i=(*n_free); i < n; i++) {
          for (j=i; j < n; j++){
            d[i] -= R[i + m*j]*x[permutation[j]];
          }
        }

        num_t lambda_perm[AS_N_U];
        uint8_t f_free = 0;
        num_t maxlam = -INFINITY;
        for (i=(*n_free); i<n; i++) {
          lambda_perm[i] = 0;
          for (j=(*n_free); j <= i; j++) {
            lambda_perm[i] += R[j + m*i]*d[j];
          }
          lambda_perm[i] *= -Ws[permutation[i]];
          if (lambda_perm[i] > maxlam) {
            maxlam = lambda_perm[i];
            f_free = i-(*n_free);
          }
        }

        if (maxlam <= AS_CONSTR_TOL) {
#ifdef AS_RECORD_COST
          if ((*iter) <= AS_RECORD_COST_N)
            costs[(*iter)-1] = calc_cost(A, b, us, n, n_v);
#endif
          exit_code = AS_SUCCESS;
          break; // feasible and optimal
        }

        // free variable
        qr_shift2(m, n, Q, R, (*n_free), (*n_free)+f_free);
        #ifdef debug_qr
        print_debug(A_ptr, Q_ptr, R_ptr, &n, &m);
        #endif

        Ws[permutation[(*n_free)+f_free]] = 0;
        uint8_t last_val = permutation[(*n_free)+f_free];
        for (i = f_free; i > 0; i--) {
          permutation[(*n_free)+i] = permutation[(*n_free)+i-1];
        }
        permutation[(*n_free)++] = last_val;

      }
    } else {

      num_t a = INFINITY;
      int i_a = 0;
      int f_bound = 0;
      int8_t i_s = 0;
      num_t temp;
      int8_t temp_s;
      for (int f=0; f < (*n_free); f++) {
        i = permutation[f];
        if (W_temp[i] == -1) {
          temp = (x[i] - lb[i]) / (x[i] - z[i]);
          temp_s = -1;
        } else if (W_temp[i] == +1) {
          temp = (ub[i] - x[i]) / (z[i] - x[i]);
          temp_s = +1;
        } else {
          continue;
        }
        if (temp < a) {
          a = temp;
          i_a = i;
          f_bound = f;
          i_s = temp_s;
        }
      }

      // update us
      for (i=0; i<n; i++) {
        num_t incr = a * (z[i] - x[i]);
        if (i == i_a) {
          x[i] = (i_s == +1) ? ub[i] : lb[i];
        } else {
          x[i] += incr;
        }
        if (x[i] != x[i]) {
        //if (isnan(us[i])) {
          // nan found
          nan_found = true;
          break;
        }
      }
      if (nan_found) {
        exit_code = AS_NAN_FOUND_US;
        break;
      }

      qr_shift2(m, n, Q, R, (*n_free)-1, f_bound);

      Ws[i_a] = i_s;
      uint8_t first_val = permutation[f_bound];
      for (i = 0; i < (*n_free)-f_bound-1; i++) {
        permutation[f_bound+i] = permutation[f_bound+i+1];
      }
      permutation[--(*n_free)] = first_val;
    }

#ifdef AS_RECORD_COST
    if ((*iter) <= AS_RECORD_COST_N)
      costs[(*iter)-1] = calc_cost(A, b, us, n, n_v);
#endif

  }
  if (exit_code == AS_ITER_LIMIT)
    (*iter)--;

  return exit_code;
}