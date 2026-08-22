#!/usr/bin/env python3
# The cusp-class assembly laws, verified against every cusp3_*_full dump:
#   (0) CONSTANCY: contrib = rho*a0 is constant across each class (every form);
#   (1) class sizes n_g = M/g;
#   (2) support: class g is DEAD iff  g == 2 mod 4  (a0 side: no eta-quotient
#       constant term)  or  N | g  (rho side: e_{eta*} entry vanishes);
#   (3) magnitude: |rho| = g/sqrt|A| = g*sqrt2/M for odd g; |rho| = g/M for 4|g;
#   (4) phase (odd g, base-independent):  8*theta_g = 1 + 6 r(g) + 4 s(g) mod 8,
#       r = #{p|g : p=3 mod 4},  s = #{p<q|g : (p/q)(q/p) = +1};
#       (4|g): rho REAL with sign mu(g/4).
#   (5) assembly: sum_g S_g = c_{eta*}(0) (checked against FORM headers).
import re, collections, cmath, math, glob, os, sys

def ff(s): return float(s.replace('\\','').replace('E','e'))

def load(path):
    txt = open(path, errors='replace').read()
    m = re.search(r'BASE (\d+) (\d+)\s+M = (\d+)', txt)
    D,N,M = int(m.group(1)), int(m.group(2)), int(m.group(3))
    rows=collections.defaultdict(list); cvals={}
    joined = txt.replace('\\\n','')
    for m2 in re.finditer(r'FORM (-?\d+):\s*\n((?:\s+c_eta.*\n)+)', txt):
        pass
    # c_eta(0)[2] line per form (eta* = component used in the dump)
    for chunk in re.split(r'FORM (-?\d+):', joined)[1:]:
        pass
    forms = re.split(r'FORM (-?\d+):', joined)
    for i in range(1, len(forms), 2):
        k = int(forms[i])
        m3 = re.search(r'c_eta\(0\)\[2\] = (\S+) \+ (\S+) i', forms[i+1])
        if m3: cvals[k] = complex(ff(m3.group(1)), ff(m3.group(2)))
    for chunk in joined.split('COSET ')[1:]:
        chunk = chunk.split('\nFORM')[0].replace('\n','')
        m2 = re.match(r'(-?\d+) c=(\d+) d=(\d+) gcd=(\d+) a0=\((\S+?),(\S+?)\) rho=\((\S+?),(\S+?)\) contrib=\((\S+?),(\S+?)\)', chunk)
        if not m2: continue
        k=int(m2.group(1)); c=int(m2.group(2)); d=int(m2.group(3)); g=int(m2.group(4))
        a0=complex(ff(m2.group(5)),ff(m2.group(6))); rho=complex(ff(m2.group(7)),ff(m2.group(8))); ct=complex(ff(m2.group(9)),ff(m2.group(10)))
        rows[(k,g)].append((c,d,a0,rho,ct))
    return D,N,M,rows,cvals

def primes(n):
    ps=[]; d=2
    while d*d<=n:
        if n%d==0:
            ps.append(d)
            while n%d==0: n//=d
        d+=1
    if n>1: ps.append(n)
    return ps

def legendre(a,p):
    a%=p
    if a==0: return 0
    return 1 if pow(a,(p-1)//2,p)==1 else -1

def mu(n):
    ps = primes(n); m = 1
    for p in ps:
        if n % (p*p) == 0: return 0
        m = -m
    return m

def theta8(g):
    ps = primes(g)
    r = sum(1 for p in ps if p%4==3)
    s = sum(1 for i in range(len(ps)) for j in range(i+1,len(ps))
            if legendre(ps[i],ps[j])*legendre(ps[j],ps[i])==1)
    return (1 + 6*r + 4*s) % 8

PASS=FAIL=0
def check(ok, msg):
    global PASS,FAIL
    if ok: PASS+=1
    else: FAIL+=1; print(f"   FAIL: {msg}")

HERE = os.path.dirname(__file__) or '.'
for path in sorted(glob.glob(os.path.join(HERE,'cusp3_*_full.out'))):
    D,N,M,rows,cvals = load(path)
    ks = sorted({k for k,_ in rows}); gs = sorted({g for _,g in rows})
    print(f"== X0^{D}({N}) M={M}: classes {gs}, forms {ks}")
    sqA = M/math.sqrt(2)
    for g in gs:
        alive_pred = not (g%4==2 or g%N==0)
        R0=None
        for k in ks:
            R = rows.get((k,g),[])
            if not R: continue
            R0 = R0 or R
            cts=[ct for *_,ct in R]
            check(max(abs(x-cts[0]) for x in cts)<1e-8, f"g={g} form {k} contrib varies")
        c,d,a0,rho,ct = sorted(R0)[0]
        check(len(R0)==M//g or None, f"g={g} n={len(R0)} != M/g={M//g}")
        alive = abs(rho)>1e-30
        check(alive==alive_pred, f"g={g}: alive={alive} but law predicts {alive_pred}")
        if not alive: continue
        if g%2==1:
            check(abs(abs(rho)-g/sqA)<1e-9, f"g={g} |rho|={abs(rho)} != g/sqrtA={g/sqA}")
            ph8 = 8*(cmath.phase(rho)/(2*math.pi)%1)
            check(abs(ph8-theta8(g))<1e-6, f"g={g} phase8={ph8:.4f} != law {theta8(g)}")
        else:
            check(abs(abs(rho)-g/M)<1e-9, f"g={g} |rho|={abs(rho)} != g/M={g/M}")
            pred = mu(g//4)*g/M
            check(abs(rho-pred)<1e-9, f"g={g} rho={rho} != mu(g/4)g/M={pred}")
    # assembly totals vs header c-values (tolerance scales with the printed
    # 10-significant-digit contribs: 33_2's classes 3/33 cancel at ~1e6)
    for k in ks:
        S = sum(sum(ct for *_,ct in rows.get((k,g),[])) for g in gs)
        mag = sum(abs(ct) for g in gs for *_,ct in rows.get((k,g),[]))
        if k in cvals:
            check(abs(S-cvals[k])<1e-6+1e-8*mag, f"form {k}: sum S_g={S} != c_eta*(0)={cvals[k]}")
print(f"\nTOTAL: {PASS} checks passed, {FAIL} failed")
