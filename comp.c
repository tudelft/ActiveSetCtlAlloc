/* 
 * Main function for testing the control allocation wls package. 
 * 
 * Till Blaha 2024
 */

#include <math.h>
#include <string.h>
#include <stdio.h>
#include "solveActiveSet.h"
#include "setupWLS.h"
#include <time.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include "test_cases.h"

static void genAbFromCtlAlloc(
    int n_v, int n_u,
    num_t G[AS_N_V*AS_N_U], num_t Wv[AS_N_V], num_t Wu[AS_N_U], num_t up[AS_N_U],
    num_t dv[AS_N_V],
    num_t A[AS_N_C*AS_N_U], num_t b[AS_N_C])
{
    num_t theta = 2.0e-9;
    num_t cond_bound = 4e5;
    num_t gamma;
    setupWLS_A(G, Wv, Wu, n_v, n_u, theta, cond_bound, A, &gamma);
    setupWLS_b(dv, up, Wv, Wu, n_v, n_u, gamma, b);
}

static void genHbetaFromAb(
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
            H[col*n + row] = 0.5 * tmp;
        }
    }

    for (row = 0; row < n; row++) {
        tmp = 0.;
        for (k = 0; k < m; k++)
            tmp -= A[row*m + k] * b[k];
        beta[row] = tmp;
    }
}

typedef enum {
    UAS,
    DAQP,
} solver_t;

#define BILLION 1000000000L

static long long unsigned int runTestCase(
    TestCase* test, num_t* us, solver_t solver, activeSetAlgoChoice algo)
{
    num_t A[AS_N_C*AS_N_U];
    num_t b[AS_N_C];
    num_t H[AS_N_U*AS_N_U];
    num_t beta[AS_N_U];
    struct timespec start, end;

    genAbFromCtlAlloc(test->n_v, test->n_u, test->JG, test->Wv, test->Wu, test->up, test->v,
                      A, b);

    clock_gettime(CLOCK_MONOTONIC, &start);
    switch(solver) {
        case UAS: {
            int Ws[AS_N_U]; memset(Ws, 0, sizeof(int)*AS_N_U);
            for (int i = 0; i < test->n_u; i++) us[i] = 0.;
            int iter;
            int n_free;
            #ifdef AS_RECORD_COST
              num_t costs[AS_RECORD_COST_N];
            #else
              num_t *costs = 0;
            #endif
            solveActiveSet(algo)(
                A, b, test->lb, test->ub, us, Ws, 100, test->n_u, test->n_v,
                &iter, &n_free, costs);
            break;
        }
        case DAQP:
            genHbetaFromAb( test->n_v+test->n_u, test->n_u, A, b, H, beta );
            break;
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    return BILLION * (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec);
}

static void verifyTestCase( TestCase* test, num_t* us, num_t* ue, num_t* ve ) {
    num_t e2 = 0;
    for (int i = 0; i < test->n_u; i++)
        e2 += (test->us[i] - us[i]) * (test->us[i] - us[i]);

    *ue = sqrt(e2 / test->n_u);

    num_t diff;
    for (int i=0; i < test->n_v; i++) {
        diff = 0.;
        for (int j=0; j < test->n_u; j++) {
            diff += test->JG[j*test->n_v + i] * (test->us[j] - us[j]);
        }
        e2 += diff*diff;
    }
    *ve = sqrt(e2 / test->n_v);
}

#define PERTURBATION 0.01
#define REPEAT_CASES 100

int main( int argc, char** argv ) {
    TestCase test_cases[N_CASES];
    fill_cases(test_cases);

    num_t ue_max = 0., ue_mean = 0.;
    num_t ve_max = 0., ve_mean = 0.;
    num_t pe_max = 0., pe_mean = 0.;

    num_t us[AS_N_U], us_perturb[AS_N_U];
    num_t ue, ve, pe;
    long long unsigned int solution_time_ns;

    for (int i = 0; i < N_CASES*REPEAT_CASES; i++) {
        int id = i % N_CASES;
        TestCase* test = &test_cases[id];

        solution_time_ns += runTestCase( test, us, UAS, AS_QR );

        if (i == id) {
            verifyTestCase( test, us, &ue, &ve);
            ue_max = (ue_max > ue) ? ue_max : ue;
            ve_max = (ve_max > ve) ? ve_max : ve;
            ue_mean += (ue - ue_mean) / (id + 1);
            ve_mean += (ve - ve_mean) / (id + 1);

            for (int j=0; j<test->n_v; j++)
                test->v[j] *= 1.f + PERTURBATION;

            runTestCase( test, us_perturb, UAS, AS_QR );

            pe = 0.;
            for (int j=0; j<test->n_u; j++) {
                pe += (us_perturb[j] - us[j]) * (us_perturb[j] - us[j]);
            }
            pe = sqrtf(pe / test->n_u);
            pe_max = (pe_max > pe) ? pe_max : pe;
            pe_mean += (pe - pe_mean) / (id + 1);
        }
    }

    printf("%llu ns mean per case\n", solution_time_ns / N_CASES / REPEAT_CASES);
    printf("ue_max: %.4f%%, ue_mean: %.4f%%\n", 100.*ue_max, 100.*ue_mean);
    printf("ve_max: %.4f%%, ve_mean: %.4f%%\n", 100.*ve_max, 100.*ve_mean);
    printf("pe_max: %.4f%%, pe_mean: %.4f%%\n", 100.*pe_max, 100.*pe_mean);

    return 0;
}
