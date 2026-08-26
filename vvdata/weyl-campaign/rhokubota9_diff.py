#!/usr/bin/env python3
"""Diff the closed-form predictions (PRED lines) against the banked 210_1
RHOV dump.  Reports per-word absolute deviation, support mismatches, and the
worst deviation."""
import gzip, sys, math

def flatten(path):
    """Join Magma's 80-column wrapped lines: a continuation line is one that
    does not start with a known keyword."""
    out = []
    for raw in open(path):
        line = raw.rstrip("\n")
        if out and not any(line.startswith(k) for k in
                           ("PRED", "DONE", "BASE", "level", "p =", "sign",
                            "Jordan", "#cosets", "x_c")):
            out[-1] += line
        else:
            out.append(line)
    return out

pred = {}       # wi -> complex or "ZERO"/"SKIP"
for line in flatten(sys.argv[1]):
    t = line.split()
    if not t or t[0] != "PRED":
        continue
    wi = int(t[1])
    if t[2] in ("ZERO", "SKIP"):
        pred[wi] = t[2]
    else:
        pred[wi] = complex(float(t[2]), float(t[3]))

meas = {}
op = gzip.open if sys.argv[2].endswith(".gz") else open
with op(sys.argv[2], "rt") as f:
    joined = []
    for raw in f:
        line = raw.rstrip("\n")
        if joined and joined[-1].startswith("RHOV") and len(joined[-1].split()) < 5 \
           and not line.startswith(("RHOV", "EMAT")):
            joined[-1] += line
        else:
            joined.append(line)
for line in joined:
    t = line.split()
    if not t or t[0] != "RHOV":
        continue
    # RHOV wi class re im
    meas[int(t[1])] = complex(float(t[3]), float(t[4]))

worst = 0.0; worst_wi = None; nok = 0; nbad = 0
support_pred_only = []; support_meas_only = []
for wi, p in sorted(pred.items()):
    m = meas.get(wi)
    if p in ("ZERO", "SKIP"):
        if m is not None and abs(m) > 1e-38:
            support_meas_only.append(wi)
        continue
    if m is None:
        if abs(p) > 1e-38:
            support_pred_only.append(wi)
        continue
    d = abs(p - m)
    if d > worst:
        worst, worst_wi = d, wi
    if d < 1e-38:
        nok += 1
    else:
        nbad += 1
        if nbad <= 12:
            print(f"BAD wi={wi} |pred-meas|={d:.3e} pred={p:.6f} meas={m:.6f} "
                  f"ratio={(m/p) if abs(p)>0 else float('nan'):.6f}")

print(f"matched exact (<1e-38): {nok}   deviating: {nbad}")
print(f"support: pred-nonzero-but-missing-in-dump {len(support_pred_only)} {support_pred_only[:8]}")
print(f"         dump-nonzero-but-pred-zero {len(support_meas_only)} {support_meas_only[:8]}")
print(f"worst dev: {worst:.3e} at wi={worst_wi}")
print(f"#pred={len(pred)} #meas={len(meas)}")
