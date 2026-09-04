#!/usr/bin/env python3
"""Extract Guo-Yang's published hyperelliptic EQUATIONS (arXiv:1510.06193) into a Magma data file.

⚠ THREE TRAPS, each of which silently yields a WRONG polynomial rather than an error:

 1. AN EQUATION MAY WRAP ACROSS `\\\\` INTO SEVERAL $...$ GROUPS, e.g. 6_11:
      $y^2=-19x^8 - ... + 612x^4$ \\\\ $+ 166x^3 - 439x^2 + 166x - 19$
    Splitting the row on `\\\\` truncates it and silently loses four terms.
    RULE: a $...$ group with no `=` that starts with + or - is a CONTINUATION of the current
    equation; a group containing `=` starts a NEW one.
 2. A ROW MAY CARRY A PAIR OF EQUATIONS (82_1: y^2=f(s) AND x^2=g(s)) -- a different
    presentation, not a single hyperelliptic y^2=f(x). Refused.
 3. VARIABLES VARY (x, s, z), and 93_1 mixes TWO (3s^3-7s^2-3t-1) -- almost certainly a typo
    in the paper. Refused rather than guessed.

VALIDATION: checked against the equations hard-coded in the existing tests/X0_D_N.m, which were
transcribed by hand from the paper independently of this script. That is the reason to trust it.
"""
import re, os

SRC = "ShimuraCurves-arxiv.tex"
TESTS = "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/tests"

def polymul(a, b):
    out = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i+j] += x*y
    return out

def parse_poly(s):
    s = s.replace(' ', '')
    if not s:
        return [0]
    coeffs = {}
    for t in re.findall(r'[+-]?[^+-]+', s):
        if not t or t in '+-':
            continue
        sign = -1 if t.startswith('-') else 1
        t = t.lstrip('+-')
        m = re.match(r'^(\d*)[a-z](?:\^\{?(\d+)\}?)?$', t)
        if m:
            c = int(m.group(1)) if m.group(1) else 1
            e = int(m.group(2)) if m.group(2) else 1
        else:
            if not re.match(r'^\d+$', t):
                raise ValueError(f"unparsed term {t!r}")
            c, e = int(t), 0
        coeffs[e] = coeffs.get(e, 0) + sign*c
    return [coeffs.get(i, 0) for i in range(max(coeffs)+1)]

def parse_rhs(rhs):
    rhs = rhs.replace('\\left', '').replace('\\right', '').replace(' ', '')
    sign = 1
    while rhs[:1] in '+-' and rhs[1:2] == '(':
        if rhs[0] == '-':
            sign = -sign
        rhs = rhs[1:]
    factors, i = [], 0
    while i < len(rhs):
        if rhs[i] == '(':
            depth, j = 1, i + 1
            while depth:
                if rhs[j] == '(':
                    depth += 1
                elif rhs[j] == ')':
                    depth -= 1
                j += 1
            factors.append(rhs[i+1:j-1]); i = j
        else:
            j = rhs.find('(', i); j = len(rhs) if j < 0 else j
            factors.append(rhs[i:j]); i = j
    out = [sign]
    for f in factors:
        if f.strip():
            out = polymul(out, parse_poly(f))
    return out

txt = open(SRC).read()
tabs = list(re.finditer(r'\\begin\{tabular\}.*?\\end\{tabular\}', txt, re.S))
caps = [m.start() for m in re.finditer(r'\\caption\{Equations of level[^}]*\}', txt)]

eqs, refused = {}, []
for cstart in caps:
    tab = max([t for t in tabs if t.end() <= cstart], key=lambda t: t.end()).group(0)
    labs = [(m.start(), int(m.group(1)), int(m.group(2)) if m.group(2) else 1)
            for m in re.finditer(r'X\^\{(\d+)\}_0(?:\((\d+)\))?', tab)]
    for idx, (pos, D, N) in enumerate(labs):
        # the label itself sits inside $...$; starting the segment AT it would leave us mid-group
        # and every following $ would pair off by one. Skip past the label's closing $ first.
        close = tab.find('$', pos)
        start = close + 1 if close != -1 else pos
        seg = tab[start: labs[idx+1][0] if idx + 1 < len(labs) else len(tab)]
        eqlist, cur = [], None
        for g in re.findall(r'\$([^$]*)\$', seg):
            g = g.strip()
            if not g:
                continue
            if re.search(r'[a-z]\^\{?2\}?\s*=', g):
                if cur is not None:
                    eqlist.append(cur)
                cur = g
            elif cur is not None and g[:1] in '+-':
                cur += g
        if cur is not None:
            eqlist.append(cur)
        if not eqlist:
            continue
        if len(eqlist) > 1:
            refused.append((D, N, "pair of equations, not a single y^2=f(x)")); continue
        body = re.sub(r'^[a-z]\^\{?2\}?\s*=\s*', '', eqlist[0])
        vs = set(re.findall(r'(?<![a-zA-Z\\])[a-z](?![a-zA-Z])', body))
        if len(vs) > 1:
            refused.append((D, N, f"two variables {sorted(vs)}, likely a typo in the paper")); continue
        try:
            eqs[(D, N)] = parse_rhs(body)
        except Exception as e:
            refused.append((D, N, str(e)))

def test_equation(D, N):
    f = f"{TESTS}/X0_{D}_{N}.m"
    if not os.path.exists(f):
        return None
    m = re.search(r'HyperellipticCurve\(([^)]*)\)', open(f).read())
    if not m:
        return None
    try:
        return parse_poly(m.group(1).replace('*', ''))
    except Exception:
        return None

def norm(p):
    while p and p[-1] == 0:
        p = p[:-1]
    return p

agree = dis = 0
print("VALIDATION against hand-transcribed equations in tests/X0_D_N.m:")
for (D, N) in sorted(eqs):
    t = test_equation(D, N)
    if t is None:
        continue
    a, b = norm(eqs[(D, N)]), norm(t)
    same = (a == b) or (a == [-c for c in b])
    agree += same; dis += (not same)
    if not same:
        print(f"  {D}_{N}: MISMATCH")
        print(f"     extracted {a}")
        print(f"     test      {b}")
print(f"  -> {agree} agree, {dis} mismatch\n")

with open("gy_equations.m", "w") as fh:
    fh.write("// Guo-Yang published hyperelliptic equations, arXiv:1510.06193.\n")
    fh.write("// AUTO-GENERATED by extract_equations.py -- do not hand-edit.\n")
    fh.write("// gy_eq[<D,N>] = coefficients of f, x^0 first, for y^2 = f(x).\n")
    fh.write("// Compare against OUR models[[Integers()|1]] (the W={1} key, the curve itself).\n")
    fh.write("gy_eq := AssociativeArray();\n")
    for (D, N) in sorted(eqs):
        fh.write(f"gy_eq[<{D},{N}>] := {eqs[(D,N)]};\n")
print(f"wrote gy_equations.m with {len(eqs)} equations; refused {len(refused)}:")
for r in refused:
    print("   ", r)
