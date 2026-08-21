from fractions import Fraction as F
import re, itertools

# required correction per form: C*S_A + 2*mult*sigma + T_f = 0, with C the global constant of the
# A-side series (common to all forms; absorb into the units of T and sigma): use C = 1 units, i.e.
# fit T_f = -S_A - 2*mult*sigma_hat  with sigma_hat = sigma/C.
SA   = {-2:F(-5,4), -1:F(-5,4), 9:F(-5), 10:F(-25,4), 11:F(-15,4), 12:F(-15,2),
        13:F(-5,4), 14:F(-15,4), 15:F(-5,4)}
mult = {-2:2, -1:4, 9:0, 10:0, 11:4, 12:2, 13:4, 14:-2, 15:2}

# positive 0-coset coefficients
pos = {}
for line in open('posdump.out'):
    m = re.match(r'POS (-?\d+) \[(.*)\]', line.strip())
    if m:
        k = int(m.group(1))
        pos[k] = [F(x.strip()) for x in m.group(2).split(',')]
ks = sorted(SA)
print(f"forms: {ks};  positive coeffs to n = {len(pos[ks[0]])}")

# Hypothesis: T_f = kappa * Sum_{r>=1} r^a * c_0(t r^2), truncated at t r^2 <= 36.
# Fit (kappa, sigma_hat) per (t, a) on two forms, test on the other seven.
best = []
for t in [1,2,3,5,6,10,15,30]:
    for a in [0,1,2]:
        # theta sums per form
        th = {k: sum(F(r)**a * pos[k][t*r*r - 1] for r in range(1,7) if t*r*r <= 36) for k in ks}
        # solve kappa*th + 2*mult*sig = -SA on forms 9, -1 (mult 0 and 4)
        k1, k2 = 9, -1
        A11, A12, b1 = th[k1], F(2*mult[k1]), -SA[k1]
        A21, A22, b2 = th[k2], F(2*mult[k2]), -SA[k2]
        det = A11*A22 - A12*A21
        if det == 0: continue
        kap = (b1*A22 - b2*A12)/det
        sig = (A11*b2 - A21*b1)/det
        resid = [k for k in ks if kap*th[k] + 2*mult[k]*sig + SA[k] != 0]
        nfit = 9 - len(resid)
        best.append((nfit, t, a, kap, sig, resid))
best.sort(reverse=True)
for nfit, t, a, kap, sig, resid in best[:8]:
    print(f"t={t:>2} a={a}: fits {nfit}/9  kappa={kap}  sigma={sig}  fails at {resid}")
