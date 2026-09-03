#!/usr/bin/env python3
"""Emit the permanent tests/_offline/GuoYang_D_N.m files for bases that passed verification."""
import json
from collections import Counter

with open("tables.json") as f:
    tables = json.load(f)

PASSING = [
    (10,11),(10,13),(10,23),(134,1),(14,3),(14,5),(146,1),(15,2),(194,1),(206,1),
    (21,2),(22,3),(22,5),(26,1),(35,1),(38,1),(39,1),(39,2),(51,1),(55,1),(57,1),
    (58,1),(6,11),(6,17),(6,19),(6,29),(6,31),(6,37),(62,1),(74,1),(82,1),(86,1),
    (87,1),(94,1),
]

TEMPLATE = '''import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^{D}({N})", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
{rows}
];
test_gy_table({D}, {N}, gy);
'''

OUTDIR = "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/tests/_offline"
count = 0
for (D, N) in PASSING:
    t = next(x for x in tables if x["D"] == D and x["N"] == N)
    disc_counts = Counter(r["disc"] for r in t["rows"])
    dup = {d for d, n in disc_counts.items() if n > 1}
    rows = []
    for r in t["rows"]:
        if r["disc"] in dup:
            continue
        v = r["vals"][0]
        if v["kind"] == "rat":
            num, den = v["val"]
            rows.append(f"<{r['disc']}, {num}>" if den == 1 else f"<{r['disc']}, {num}/{den}>")
        elif v["kind"] == "oo":
            rows.append(f"<{r['disc']}, Infinity()>")
    content = TEMPLATE.format(D=D, N=N, rows=",\n".join(rows))
    fname = f"{OUTDIR}/GuoYang_{D}_{N}.m"
    with open(fname, "w") as f:
        f.write(content)
    count += 1

print(f"wrote {count} test files")
