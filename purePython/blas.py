from math import sqrt
from __main__ import CYCLES

def daxpy( n, a, x, y ):
    global CYCLES
    for i in range(n):
        y[i] += a*x[i]
        CYCLES += 1

def dscal( n, a, x ):
    global CYCLES
    for i in range(n):
        x[i] *= a
        CYCLES += 1

def ddot( n, x, y ):
    global CYCLES
    out = x[0]*y[0]
    CYCLES += 1
    for i in range(1, n):
        out += x[i]*y[i]
        CYCLES += 1
    return out

def dnrm2( n, x ):
    global CYCLES
    if n <= 0:
        return 0
    out = x[0]*x[0]
    CYCLES += 1
    for i in range(1, n):
        CYCLES += 1
        out += x[i]*x[i]
    CYCLES += 14
    return sqrt(out)
