# Local Whittaker at the LEVEL PRIME N, by the alpha_k / G(X) recipe the m=0 paper uses:
#   alpha_k = p^{-k(n-1)} #{ x in (mu+L)/p^k L : Q(x) = m mod p^k },  n = 2
#   W_{m,p}(s,phi) = (1-X) G(X),  G(X) = sum_k alpha_k X^k,  X = p^{-s}
#
# At p = N the lattice is the level-N hyperbolic plane: L = Z e + Z f, Q(ue+vf) = N u v.
# Then L^v/L = (Z/N)^2 with Q(a/N, b/N) = ab/N mod 1, so isotropic <=> ab = 0 mod N,
# giving exactly 2N-1 isotropic cosets -- which is MEASURED FACT 1 (2N-1, 12 bases).
from fractions import Fraction

def alpha(N, k, m, mu):          # mu = (a,b) in (Z/N)^2 representing (a/N, b/N)
    a, b = mu; Nk = N**k
    cnt = 0
    for s in range(Nk):
        for t in range(Nk):
            # x = (a/N + s) e + (b/N + t) f ; Q(x) = N(a/N+s)(b/N+t) = ab/N + a t + b s + N s t
            # work N*Q(x) to stay integral, compare mod N^{k+1}
            NQ = a*b + N*(a*t + b*s + N*s*t)
            if (NQ - N*m) % (N*Nk) == 0: cnt += 1
    return Fraction(cnt, N**k)   # n-1 = 1

def Wcoeffs(N, m, mu, K=4):
    al = [alpha(N, k, m, mu) for k in range(K)]
    return [al[0]] + [al[k]-al[k-1] for k in range(1, K)]   # (1-X)G(X)

for N in (2,3,5):
    print(f"--- N={N} ---")
    print(f"  mu=0,   m=0 : W coeffs {Wcoeffs(N,0,(0,0))}   (paper: 1,(N-1),(N-1),... = 1,{N-1},{N-1},...)")
    for mu in [(1,0),(0,1)] + ([(2,0)] if N>2 else []):
        for m in (0,1,2,3):
            print(f"  mu={mu}, m={m} : W coeffs {Wcoeffs(N,m,mu)}")

print()
print("=== ZERO coset, varying m: does the level prime see ord_N(m)? ===")
for N in (2,3,5):
    for m in range(0,7):
        c = Wcoeffs(N,m,(0,0),K=4)
        print(f"  N={N} mu=0 m={m} (ord_N={0 if m==0 else next(k for k in range(9) if m%N**(k+1))}) W = {[str(x) for x in c]}")
