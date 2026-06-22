#!/usr/bin/env python3
"""Surgically rewrite the TestInWhichProved string of each curve record in
data/curves_after_UpdateCurves8.dat to the corrected attribution from
data/attribution_map.txt (keyed by CurveID), leaving every other byte untouched.
This produces a minimal git diff: only the 752 mis-attributed lines change."""
import re, sys

DAT = "data/curves_after_UpdateCurves8.dat"
MAP = "data/attribution_map.txt"

mapping = {}
entry_re = re.compile(r'^(\d+)\t(.*)$')
cur = None
with open(MAP) as f:
    for line in f:
        line = line.rstrip("\n")
        m = entry_re.match(line)
        if m:                              # new "id<TAB>value" entry
            cur = int(m.group(1))
            mapping[cur] = m.group(2)
        elif line and cur is not None:     # wrapped continuation of previous value
            mapping[cur] += line

with open(DAT) as f:
    text = f.read()

# Long records are line-wrapped in the .dat, so allow whitespace (incl. newlines)
# wherever Magma may have broken the token.
token = re.compile(r'<"CurveID",\s*(\d+)>|<"TestInWhichProved",\s*"([^"]*)">')

cur_id = None
out = []
last = 0
changes = 0
for m in token.finditer(text):
    if m.group(1) is not None:           # CurveID token
        cur_id = int(m.group(1))
        continue
    # TestInWhichProved token
    old = m.group(2)
    if cur_id is None or cur_id not in mapping:
        sys.exit(f"no mapping for CurveID {cur_id} near offset {m.start()}")
    new = mapping[cur_id]
    if new != old:
        out.append(text[last:m.start()])
        out.append(f'<"TestInWhichProved", "{new}">')
        last = m.end()
        changes += 1
out.append(text[last:])
result = "".join(out)

with open(DAT, "w") as f:
    f.write(result)
print(f"rewrote {changes} TestInWhichProved values in {DAT}")
