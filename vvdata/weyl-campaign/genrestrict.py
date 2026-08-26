#!/usr/bin/env python3
# genrestrict.py <mono dump> <eisrat log> <pairs> <out .m>
#   pairs: "7,1;7,30"  -- restrict the Gross-span fit to EXACTLY these (D',R),
#   i.e. the gauge check that a competing support does or does not carry E_eis.
# Generalizes the hand-written restrict210A/B.m to any base.
import re, sys, os
mono, eislog, pairstr, outp = sys.argv[1:5]
mt = open(mono, errors='replace').read()
mb = re.search(r'BASE (\d+) (\d+)\s+M = (\d+)\s+ds = \[(.*?)\]', mt, re.S)
D, N = int(mb.group(1)), int(mb.group(2))
ds = [int(x) for x in mb.group(4).replace('\n', '').split(',')]
DN = D * N
et = open(eislog, errors='replace').read()
assert 'RATIONAL E_eis FOUND' in et, "eisrat did not certify a rational E_eis"
rs, cs = [], []
for m in re.finditer(r'beta\[\d+\] = (\S+)\s+r = \[(.*?)\]', et, re.S):
    c = m.group(1)
    if c == '0':
        continue
    cs.append(c)
    rs.append([int(x) for x in m.group(2).replace('\n', '').split(',')])
assert rs, "no nonzero beta rows"

here = os.path.dirname(os.path.abspath(sys.argv[0]))
tmpl = None
for cand in (os.path.join(here, 'allgross_template.m'),
             '/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/'
             'dbce22ad-34ec-4786-94b9-f6da57354d63/scratchpad/sci/allgross_template.m'):
    if os.path.exists(cand):
        tmpl = open(cand).read()
        break
assert tmpl, "allgross_template.m not found"

pairs = [tuple(int(x) for x in p.split(',')) for p in pairstr.split(';')]
plit = '[ ' + ', '.join(f'<{a}, {b}>' for a, b in pairs) + ' ]'
# replace the template's automatic pair enumeration (a NESTED for, so match
# through to the line that consumes `pairs`) with the fixed list
tmpl = re.sub(r'pairs := \[\];\n[\s\S]*?(?=printf "definite structures)',
              f'pairs := {plit};\n', tmpl, count=1)
assert 'pairs := [ <' in tmpl, "pair-list substitution failed"

rows = ',\n        '.join('[' + ','.join(str(x) for x in r) + ']' for r in rs)
body = (tmpl.replace('@DN@', str(DN)).replace('@DSV@', str(ds).replace(' ', ''))
            .replace('@RS@', rows).replace('@CS@', ', '.join(cs))
            .replace('@TAG@', f"{D}_{N}"))
open(outp, 'w').write(body)
print(f"# wrote {outp}: {len(rs)} monomials, DN={DN}, restricted to {plit}")
