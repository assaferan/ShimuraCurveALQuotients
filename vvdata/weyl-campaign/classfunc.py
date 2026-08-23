#!/usr/bin/env python3
# The per-class functional: for each cusp class g on each cusp3-dumped base,
# solve S_g(f) = sum_m W_g(m) c_f(-m) over the oo principal-part channels
# (rref over the 9-form panel), report forced combinations, and test the
# Ramanujan-sum ansatz W_g(m) = alpha_g * c_q(m) for q in divisors(M).
#
# Data: cusp3_*_full.out (class sums, 10 sig digits) + pp3/pp2 dumps (exact pp).
import re, collections, cmath, math, glob, os, sys
from fractions import Fraction

CAMP = sys.argv[1] if len(sys.argv) > 1 else '.'

def ff(s): return float(s.replace('\\','').replace('E','e'))

def load_cusp3(path):
    txt = open(path, errors='replace').read()
    m = re.search(r'BASE (\d+) (\d+)\s+M = (\d+)', txt)
    D,N,M = int(m.group(1)), int(m.group(2)), int(m.group(3))
    rows = collections.defaultdict(list)
    joined = txt.replace('\\\n','')
    for chunk in joined.split('COSET ')[1:]:
        chunk = chunk.split('\nFORM')[0].replace('\n','')
        m2 = re.match(r'(-?\d+) c=(\d+) d=(\d+) gcd=(\d+) a0=\((\S+?),(\S+?)\) rho=\((\S+?),(\S+?)\) contrib=\((\S+?),(\S+?)\)', chunk)
        if not m2: continue
        k=int(m2.group(1)); g=int(m2.group(4))
        ct=complex(ff(m2.group(9)),ff(m2.group(10)))
        rows[(k,g)].append(ct)
    S = {}
    for (k,g),cts in rows.items():
        S[(k,g)] = sum(cts)
    return D,N,M,S

def load_pp(D,N):
    # pp3 preferred, else pp2; channels: ('oo', m) exact Fraction coeff;
    # 0-side channels ('0', m) kept separately (weights allowed but reported).
    for stem in (f'pp3_{D}_{N}.out', f'pp2_{D}_{N}.out'):
        p = os.path.join(CAMP, stem)
        if os.path.exists(p): break
    else:
        return None
    pp = collections.defaultdict(dict)
    txt = open(p, errors='replace').read().replace('\\\n','')
    for line in txt.splitlines():
        t = line.split()
        if not t or t[0] != 'COEF': continue
        _, d_, n_, k, cusp, m1, m2, c = t
        if int(d_) != D or int(n_) != N: continue
        pp[int(k)][(cusp, m1)] = Fraction(c)
    return pp

def primes(n):
    ps=[]; d=2
    while d*d<=n:
        if n%d==0:
            ps.append(d)
            while n%d==0: n//=d
        d+=1
    if n>1: ps.append(n)
    return ps

def mu(n):
    m=1
    for p in primes(n):
        if n%(p*p)==0: return 0
        m=-m
    return m

def phi(n):
    r=n
    for p in primes(n): r = r//p*(p-1)
    return r

def ramanujan(q, m):
    # c_q(m) = mu(q/gcd)*phi(q)/phi(q/gcd)
    if q == 1: return 1
    d = math.gcd(q, m)
    qd = q//d
    mq = mu(qd)
    if mq == 0: return 0
    return mq*phi(q)//phi(qd)

def divisors(n):
    ds = [d for d in range(1,n+1) if n%d==0]
    return ds

def lstsq_1d(xs, ys):
    # fit ys ~ alpha*xs (complex ys, real xs); return alpha, max residual
    num = sum(y*x for x,y in zip(xs,ys)); den = sum(x*x for x in xs)
    if den == 0:
        alpha = 0
    else:
        alpha = num/den
    res = max(abs(y-alpha*x) for x,y in zip(xs,ys))
    return alpha, res

for path in sorted(glob.glob(os.path.join(CAMP,'cusp3_*_full.out'))):
    D,N,M,S = load_cusp3(path)
    pp = load_pp(D,N)
    if pp is None:
        print(f"== X0^{D}({N}): NO pp dump, skipped"); continue
    ks = sorted({k for k,_ in S}); gs = sorted({g for _,g in S})
    scale = max(abs(v) for v in S.values())
    print(f"\n== X0^{D}({N}) M={M} classes={gs} forms={ks} (S scale {scale:.2e})")
    # channel lists
    oo_ms = sorted({int(m) for k in ks for (cusp,m) in pp.get(k,{}) if cusp=='oo'})
    zero_ch = sorted({m for k in ks for (cusp,m) in pp.get(k,{}) if cusp=='0'})
    if zero_ch: print(f"   0-side channels present: {zero_ch}")
    for g in gs:
        ys = [S.get((k,g), 0) for k in ks]
        if max(abs(y) for y in ys) < 1e-6*scale:  # dead class
            continue
        imax = max(abs(y.imag) for y in ys)
        tag = "" if imax < 1e-6*scale else f"  [IMAG {imax:.2e}]"
        # Ramanujan ansatz over oo channels, all moduli q | M
        best = []
        for q in divisors(M):
            xs = [float(sum(pp[k].get(('oo',str(m)),0)*ramanujan(q,m) for m in oo_ms)) for k in ks]
            alpha, res = lstsq_1d(xs, [y.real for y in ys])
            if res < 1e-5*max(scale,1): best.append((q, alpha, res))
        line = f"  g={g:>3} S_g/forms: " + " ".join(f"{y.real:+.4g}" for y in ys) + tag
        print(line)
        if best:
            for q,alpha,res in best:
                print(f"        ANSATZ HIT: W_g(m) = {alpha:+.6g} * c_{q}(m)   (max resid {res:.1e})")
        else:
            print(f"        no single-Ramanujan fit (moduli q|{M})")
print("\ndone")
