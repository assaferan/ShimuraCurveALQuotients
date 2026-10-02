#!/usr/bin/env python3
"""Fit the Weil-stage cost against the measured times.

    python3 fit.py retime.psv weil_design.log

retime.psv is one line per finished curve, `id|D,N,g,#W,Qmax|total cpu s|primes done`, scraped from
logs/weil_main2_<id>.log; weil_design.log lists each curve's admissible primes (weil_design.m).
Compares three candidate estimates and fits power laws.  Conclusion, 2026-10-02 on 34 curves:
#W * sum p^(g/2) orders 216 of the 220 pairs that differ by more than 10x in time, the sum of p^g
only 166; the level does not enter (exponent -0.07); the exact term count with the Q_w^(-1/2)
weights of the trace formula is WORSE than the unweighted #W version, so the cost grows with #W
rather than shrinking with the Q_w.
"""
import sys, re, math

def load(psv, design):
    times = {}
    for line in open(psv):
        i, shape, tot, np_ = line.strip().split("|")
        D, N, g, W, Q = [int(x) for x in shape.split(",")]
        times[int(i)] = dict(D=D, N=N, g=g, W=W, Q=Q, t=float(tot), np=int(np_))
    rows, cur = [], ""
    for line in open(design):
        line = line.rstrip()
        if re.match(r"^\d+\s", line):
            if cur:
                rows.append(cur)
            cur = line
        elif cur and not line.startswith(("id", "63 ")):
            cur += " " + line.strip()
    if cur:
        rows.append(cur)
    prim = {}
    for r in rows:
        m = re.match(r"(\d+)\s+(\S+)\s+(\d)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+\[(.*?)\]", r)
        if m:
            prim[int(m[1])] = [int(x) for x in m[8].split(",") if x.strip()]
    out = []
    for i, v in times.items():
        ps = prim.get(i)
        if not ps:
            continue
        ps = ps[:v["np"]]
        out.append((i, v, v["W"] * sum(p ** (v["g"] / 2) for p in ps), sum(p ** v["g"] for p in ps)))
    return out

def inversions(data, key, thr):
    bad = [(a[0], b[0]) for a in data for b in data if a[1]["t"] > thr * b[1]["t"] and a[key] < b[key]]
    tot = sum(1 for a in data for b in data if a[1]["t"] > thr * b[1]["t"])
    return len(bad), tot, bad

def fit(data, cols):
    A = [[1.0] + [math.log(f(v, s2)) for f in cols] for i, v, s2, s1 in data]
    y = [math.log(v["t"]) for i, v, s2, s1 in data]
    n = len(cols) + 1
    M = [[sum(A[k][r] * A[k][s] for k in range(len(A))) for s in range(n)] for r in range(n)]
    b = [sum(A[k][r] * y[k] for k in range(len(A))) for r in range(n)]
    for r in range(n):
        p = max(range(r, n), key=lambda q: abs(M[q][r]))
        M[r], M[p] = M[p], M[r]
        b[r], b[p] = b[p], b[r]
        for q in range(r + 1, n):
            f = M[q][r] / M[r][r]
            for s in range(r, n):
                M[q][s] -= f * M[r][s]
            b[q] -= f * b[r]
    x = [0] * n
    for r in reversed(range(n)):
        x[r] = (b[r] - sum(M[r][s] * x[s] for s in range(r + 1, n))) / M[r][r]
    pred = [math.exp(sum(A[k][r] * x[r] for r in range(n))) for k in range(len(A))]
    rat = [math.exp(y[k]) / pred[k] for k in range(len(A))]
    return x, max(rat) / min(rat)

if __name__ == "__main__":
    data = load(sys.argv[1], sys.argv[2])
    print(f"{len(data)} curves with a measured time\n")
    for thr in (2, 3, 5, 10):
        for name, k in (("#W sum p^(g/2)", 2), ("sum p^g", 3)):
            n, tot, bad = inversions(data, k, thr)
            print(f"  times differing >{thr:>2}x: {name:<16} mis-orders {n:>3} of {tot:>3}"
                  + (f"   {bad}" if bad and len(bad) <= 6 else ""))
        print()
    for label, cols in (
        ("#W^a (sum p^(g/2))^b", [lambda v, s: v["W"], lambda v, s: s]),
        ("... * (D N)^c", [lambda v, s: v["W"], lambda v, s: s, lambda v, s: v["D"] * v["N"]]),
    ):
        x, spread = fit(data, cols)
        print(f"  {label:<22} exponents {[round(e, 2) for e in x[1:]]}  residual spread {spread:.0f}x")
