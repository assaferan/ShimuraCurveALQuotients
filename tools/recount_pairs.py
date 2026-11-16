"""Recount the pairs (D, N), D > 1, that FindPairs keeps, written independently of the Magma code.

A pair is kept unless phi(D) psi(N) / 2^omega(DN) > 24 (p^2 + 1)/(p - 1) for p the smallest prime not
dividing DN: at a prime of good reduction the supersingular points number at least
phi(D)(p-1)psi(N)/12 [Ribet 1989, Thm 3.4], and a curve of gonality 2 has at most 2(p^2 + 1) points
over F_{p^2}.  The search range N <= C/phi(D), C = 2^6 * 435, is the one FindPairs uses.
Expected output: 2717 pairs, maximal D 99330, 505 values of D.
"""
from fractions import Fraction
from math import gcd, ceil
N_MAX=300000
spf=list(range(N_MAX+1))
for i in range(2,int(N_MAX**.5)+1):
    if spf[i]==i:
        for j in range(i*i,N_MAX+1,i):
            if spf[j]==j: spf[j]=i
def fac(n):
    f={}
    while n>1:
        p=spf[n]; f[p]=f.get(p,0)+1; n//=p
    return f
def phi(n):
    r=n
    for p in fac(n): r=r//p*(p-1)
    return r
def psi(n):
    r=Fraction(n)
    for p in fac(n): r*=Fraction(p+1,p)
    return r
PRIMES=[q for q in range(2,200) if all(q%r for r in range(2,q))]
g=lambda p: Fraction(24*(p*p+1),p-1)
C=64*435
pairs=[]
for D in range(1,N_MAX+1):
    f=fac(D)
    if any(e>1 for e in f.values()) or len(f)%2: continue
    lD=Fraction(phi(D),2**len(f))
    Nmax=ceil(Fraction(C,phi(D)))
    for N in range(1,Nmax+1):
        if gcd(D,N)!=1: continue
        p=next(q for q in PRIMES if (D*N)%q)
        fN=fac(N); lN=psi(N)/2**len(fN)
        if D>1:
            if lD*lN<=g(p): pairs.append((D,N))
        else:
            # D = 1: add the cusps as the code does (LowerBound); reproduce only the D > 1 count here
            pass
print("D>1 pairs:",len(pairs),"max D",max(d for d,n in pairs),"values of D (D>1)",len(set(d for d,n in pairs)))
