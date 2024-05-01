CYCLES = 0
EPS = 1e-15

from blas import dnrm2, daxpy, dscal, ddot
from math import hypot

def dorg2r( m, n, k, U, tau ):
    global CYCLES
    for ik in range(n):
    #    A[:k, k] = 0.
        U[ik, ik] = 1.
        CYCLES += 1

    # fully optimized
    ik = k-1
    for iin in range(n-1,k-2,-1):
        tmp = U[iin, k-1]
        for im in range(k-1, m):
            U[im, iin] = -tau[k-1] * U[im, k-1] * tmp
            if im == iin:
                U[im, iin] += 1.
            CYCLES += 2

    ik -= 1
    while ik >= 0:
        iin = n - 1
        while iin >= ik:
            tmp = 0
            if iin == ik:
                tmp = 1
            else:
                #print(A[k+1:, k], A[k+1:, m])
                for im in range(ik+1, m):
                #for j in range(k+1, k+nmb+1):
                    tmp += U[im, ik] * U[im, iin]
                    CYCLES += 1
            tmp *= tau[ik]
            CYCLES += 1
            for im in range(ik, m):
            #for j in range(k, k+nmb+1):
                CYCLES += 1
                if ik == iin:
                    U[im, iin] = -U[im, ik] * tmp
                    if im == iin:
                        U[im, iin] += 1.
                else:
                    #print(A[k:, m], A[k:, k], tmp)
                    if im == ik:
                        U[im, iin] = -tmp
                    else:
                        U[im, iin] -= U[im, ik] * tmp
            iin -= 1
        ik -= 1

def dorg2r_sparse( m, n, mb, k, U, tau ):
    global CYCLES
    for ik in range(n):
    #    A[:k, k] = 0.
        U[ik, ik] = 1.
        CYCLES += 1

    # fully optimized
    ik = k-1
    for iin in range(n-1,k-2,-1):
        tmp = U[iin, k-1]
        for im in range(k-1, m):
            U[im, iin] = -tau[k-1] * U[im, k-1] * tmp
            if im == iin:
                U[im, iin] += 1.
            CYCLES += 2

    ik -= 1
    while ik >= 0:
        iin = n - 1
        while iin >= ik:
            tmp = 0
            if iin == ik:
                tmp = 1
            else:
                for im in range(ik+1, ik+mb+1):
                    tmp += U[im, ik] * U[im, iin]
                    CYCLES += 1
            tmp *= tau[ik]
            CYCLES += 1
            for im in range(ik, ik+mb+1):
                CYCLES += 1
                if ik == iin:
                    U[im, iin] = -U[im, ik] * tmp
                    if im == iin:
                        U[im, iin] += 1.
                else:
                    if im == ik:
                        U[im, iin] = -tmp
                    else:
                        U[im, iin] -= U[im, ik] * tmp
            iin -= 1
        ik -= 1

def dlarfg( n, x ):
    global CYCLES, EPS, SQRT2
    nu = dnrm2( n-1, x[1:] ); CYCLES += (14 + n-1) if n-1 > 0. else 0.

    if nu > EPS:
        #https://www.netlib.org/lapack/explore-html/d8/d0d/group__larfg_gadc154fac2a92ae4c7405169a9d1f5ae9.html#gadc154fac2a92ae4c7405169a9d1f5ae9
        beta = ( -1. if x[0] >= 0 else +1. ) * hypot( x[0], nu ); CYCLES += 15
        tau = (beta - x[0]) / beta;  CYCLES += 14
        inu = 1. / (x[0] - beta);  CYCLES += 14
        x[0] = beta

        dscal( n-1, inu, x[1:] ); CYCLES += n-1
        return tau
    else:
        return 0

def dlarf_left( m, n, v, tau, c ):
    global CYCLES
    if tau == 0.:
        return
    for col in range(1,n):
        tmp = tau * ddot( m, v, c[0:m, col] );  CYCLES += 1 + m
        for row in range(m):
            # potentially, col=0 is trivial and doesnt need this loop
            c[row, col] -= v[row] * tmp;  CYCLES += 1

#from math import abs

def dgeqr2( m, n, a, tau ):
    global CYCLES
    for i in range(n):
        tau[i] = dlarfg( m-i, a[i:m, i] )
        aii = a[i, i]
        a[i, i] = 1.
        dlarf_left( m-i, n-i, a[i:m, i], tau[i], a[i:m, i:n] )
        a[i, i] = aii

def dgeqr2_sparse( m, n, mb, a, tau ):
    global CYCLES
    for i in range(n):
        tau[i] = dlarfg( mb+1, a[i:m, i] )
        aii = a[i, i]
        a[i, i] = 1.
        dlarf_left( mb+1, n-i, a[i:m, i], tau[i], a[i:m, i:n] )
        #dlarf_left( m-i, n-i, a[i:m, i], tau[i], a[i:m, i:n] )
        a[i, i] = aii

if __name__=="__main__":
    import numpy as np
    np.random.seed(1)
    mb = 4
    n = 6
    m = n + mb
    A = np.triu(np.random.random((m,n)) - 0.5, -mb)
    #m = 4
    #n = 3
    #A = np.random.random((m,n)) - 0.5

    U, TAU = np.linalg.qr(A, 'raw')
    U = U.T
    Q, R = np.linalg.qr(A, 'complete')

    Awork = A.copy()
    tau = np.empty(n)
    work = np.empty(m)
    dgeqr2_sparse( m, n, mb, Awork, tau )
    #dgeqr2( m, n, Awork, tau )


    print("Sparse householder: ", CYCLES)
    CYCLES_Householder = CYCLES
    print(np.max(np.abs(U - Awork)))

    Rwork = np.triu(Awork)
    Qwork = np.zeros((m, m))
    Qwork[:, :n] = Awork
    dorg2r_sparse( m, m, mb, np.sum(tau > 0.), Qwork, tau )
    #dorg2r( m, m, np.sum(tau > 0.), Qwork, tau )

    print("Sparse org2r: ", CYCLES - CYCLES_Householder)
    print("Total: ", CYCLES)
    print(np.max(np.abs(Qwork - Q)))

    Atest = np.array([
        [ -80, -36.9951553, 64.0774918, 0.0500000007, 0.000101227481, 0  , 0  , 0  , 0  , 0  , ],
        [ -80, -82.7360229, 6.32455531e-15, -0.0500000007, 0  , 0.00015372863, 0  , 0  , 0  , 0  , ],
        [ -80, -39.7080421, -68.7763443, 0.0500000007, 0  , 0  , 0.000136507821, 0  , 0  , 0  , ],
        [ -80, 38.6960945, -67.0236053, -0.0500000007, 0  , 0  , 0  , 0.000143456768, 0  , 0  , ],
        [ -80, 70.6733932, 5.05964425e-14, 0.0500000007, 0, 0, 0, 0, 9.8791199e-05, 0, ],
        [ -80, 35.3365135, 61.2046356, -0.0500000007, 0, 0, 0, 0, 0, 0.000160468131, ],
    ]).T

    dgeqr2( 10, 6, Atest, tau )

