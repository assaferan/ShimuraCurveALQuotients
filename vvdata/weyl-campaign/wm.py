# Exact q-expansion of E_eis on 15_2 and the functional W(m) = -a(m).
from fractions import Fraction
DEPTH = 45
def mulpoly(a, b):
    c = [Fraction(0)]*DEPTH
    for i, x in enumerate(a):
        if x == 0: continue
        for j, y in enumerate(b):
            if i+j >= DEPTH: break
            if y != 0: c[i+j] += x*y
    return c
def eta_unit_pow(d, e):
    ser = [Fraction(0)]*DEPTH; ser[0] = Fraction(1)
    n = 1
    while d*n < DEPTH:
        if e > 0:
            fac = [Fraction(0)]*DEPTH; fac[0] = Fraction(1)
            fac[d*n] = Fraction(-1)
            for _ in range(e): ser = mulpoly(ser, fac)
        else:
            inv = [Fraction(0)]*DEPTH
            k = 0
            while d*n*k < DEPTH: inv[d*n*k] = Fraction(1); k += 1
            for _ in range(-e): ser = mulpoly(ser, inv)
        n += 1
    return ser
def mono(r, ds):
    s1 = sum(d*e for d, e in zip(ds, r))
    assert s1 % 24 == 0
    shift = s1 // 24
    ser = [Fraction(0)]*DEPTH; ser[0] = Fraction(1)
    for d, e in zip(ds, r):
        if e: ser = mulpoly(ser, eta_unit_pow(d, e))
    out = [Fraction(0)]*DEPTH
    for i, x in enumerate(ser):
        if i+shift < DEPTH: out[i+shift] = x
    return out
ds = [1,2,3,4,5,6,10,12,15,20,30,60]
P1 = mono([-6,15,0,-6,0,0,0,0,0,0,0,0], ds)
P3 = mono([-3,7,0,-3,1,0,0,0,0,1,0,0], ds)
P5 = mono([-1,0,0,4,-1,0,3,0,0,-2,0,0], ds)
P7 = mono([-1,2,0,0,0,0,0,0,-3,0,7,-2], ds)
E = [Fraction(4,5)*a - 4*b - Fraction(4,5)*c + 4*d for a,b,c,d in zip(P1,P3,P5,P7)]
print("m : a_E(m)")
for m in range(0, DEPTH):
    print(f"{m:3} : {E[m]}")
