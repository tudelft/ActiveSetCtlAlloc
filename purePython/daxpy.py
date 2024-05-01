

def daxpy( n, a, x, y ):
    for i in range(n):
        y[i] += a*x[i]