
def dorg2r( m, n, k, A, TAU ):
    # fully optimized
    ik = k-1
    for im in range(m-1,k-2,-1):
        tmp = A[im, k-1]
        for j in range(k-1, m):
            A[j, im] = -TAU[k-1] * A[j, k-1] * tmp
            if j == im:
                A[j, im] += 1.
            MAC_COUNT += 2

    #print(Q[(nk-1):, nk-1])
    #Q[(nk-1):, nk-1] -= TAU[nk-1] * A[(nk-1):, nk-1]
    #MAC_COUNT += (nm-nk)
    ik -= 1

    while ik >= 0:
        im = m - 1
        while im >= ik:
            tmp = 0
            if im == ik:
                tmp = 1
            else:
                #print(A[k+1:, k], A[k+1:, m])
                for j in range(ik+1, m):
                #for j in range(k+1, k+nmb+1):
                    tmp += A[j, ik] * A[j, im]
                    MAC_COUNT += 1
            tmp *= TAU[ik]
            MAC_COUNT += 1
            for j in range(ik, m):
            #for j in range(k, k+nmb+1):
                MAC_COUNT += 1
                if ik == im:
                    A[j, im] = -A[j, ik] * tmp
                    if j == im:
                        A[j, im] += 1.
                else:
                    #print(A[k:, m], A[k:, k], tmp)
                    if j == ik:
                        A[j, im] = -tmp
                    else:
                        A[j, im] -= A[j, ik] * tmp
            im -= 1
        ik -= 1

