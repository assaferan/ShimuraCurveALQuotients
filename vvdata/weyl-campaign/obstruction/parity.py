#!/usr/bin/env python3
# Tabulate the OBSTPAIR lines of an instrumented Borcherds search: per failing key, the obstruction
# dimension, the parity pattern of the pairings over the anchor triples, and whether some triple has
# every pairing even (the even-divisor escape hatch applies exactly then).
import re, sys, os
from collections import defaultdict, Counter
here = os.path.dirname(os.path.abspath(__file__))
pat = re.compile(r"OBSTPAIR key (\S+) infty (\S+) others \[(.*?)\] dim (\d+) pairing (\S+) parity (-?\d+) phi (.*)$", re.M)
for b in sys.argv[1:]:
    txt = open(os.path.join(here, f"obstpair_{b}.log")).read()
    trip = defaultdict(list); dims = {}; phis = defaultdict(set); pairings = defaultdict(set)
    for m in pat.finditer(txt):
        key, inf, oth, dim, pr, par, phi = m.groups()
        trip[(key, inf, oth)].append(int(par)); dims[key] = int(dim); phis[key].add(phi[:80]); pairings[key].add(pr)
    print(f"== {b}: failing keys {sorted(dims, key=int)}, obstruction dimension per key {dims}")
    for key in sorted(dims, key=int):
        ts = {t: v for t, v in trip.items() if t[0] == key}
        alleven = [t for t, v in ts.items() if all(p == 0 for p in v)]
        nonint = [t for t, v in ts.items() if any(p == -1 for p in v)]
        pc = Counter(tuple(sorted(v)) for v in ts.values())
        print(f"  key {key}: {len(ts)} failing triples; parity patterns {dict(pc)}; ALL-EVEN triples {len(alleven)}; non-integral pairings {len(nonint)}; distinct pairing values {len(pairings[key])}; distinct phi {len(phis[key])}")
        if alleven:
            print("    first all-even triple:", alleven[0])
        print("    pairing values:", sorted(pairings[key], key=lambda s: abs(float(s.replace('/', '/')) if '/' not in s else abs(eval(s))))[:20])
    res = re.search(r"OBSTPAIR RESULT.*", txt)
    print("  ", res.group(0)[:110] if res else "no result line")
    for key in sorted(dims, key=int):
        print(f"  phi at key {key}:", sorted(phis[key])[0][:80])
