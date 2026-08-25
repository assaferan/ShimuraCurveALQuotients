#!/usr/bin/env python3
# Generate a SPANNING-SUBSET polymake_solution file for rung n+k from a complete
# rung-n file plus weight-0 shift vectors (validated at M=308: rank 108/108).
#   tshift_gen.py <sol_n file> <weight0 file> <out file>
# Every emitted vector is a genuine lattice point of the larger polytope
# (poles add at oo; holomorphy elsewhere is preserved by addition).
import re, sys
def load(p):
    return [ [int(x) for x in m.split(',')] for m in re.findall(r'\\?\[ ([-\d, ]+) \]', open(p).read()) ]
sol, w0, outp = load(sys.argv[1]), load(sys.argv[2]), sys.argv[3]
shifts = [ t for t in w0 if any(t) ]
pts = { tuple(r) for r in sol }
for t in shifts:
    for r in sol:
        pts.add(tuple(a + b for a, b in zip(r, t)))
pts = sorted(pts)
print(f"# {len(sol)} base points, {len(shifts)} shifts -> {len(pts)} candidates", file=sys.stderr)
with open(outp, 'w') as f:
    f.write('[ PowerSequence(IntegerRing()) |\n')
    f.write(',\n'.join('[ ' + ', '.join(str(x) for x in r) + ' ]' for r in pts))
    f.write('\n]\n')
