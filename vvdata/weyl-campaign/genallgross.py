#!/usr/bin/env python3
# genallgross.py <mono dump> <eisrat log> <out .m>  -- build the Gross-span driver
import re, sys
mono, eislog, outp = sys.argv[1], sys.argv[2], sys.argv[3]
mt = open(mono, errors='replace').read()
mb = re.search(r'BASE (\d+) (\d+)\s+M = (\d+)\s+ds = \[(.*?)\]', mt)
D, N, M = int(mb.group(1)), int(mb.group(2)), int(mb.group(3))
ds = [int(x) for x in mb.group(4).split(',')]
DN = D*N
et = open(eislog, errors='replace').read()
assert 'RATIONAL E_eis FOUND' in et, "eisrat did not certify a rational E_eis"
rs, cs = [], []
for m in re.finditer(r'beta\[\d+\] = (\S+)\s+r = \[(.*?)\]', et):
    c = m.group(1)
    if c == '0': continue
    cs.append(c)
    rs.append([int(x) for x in m.group(2).split(',')])
assert rs, "no nonzero beta rows"
rows = ',\n        '.join('[' + ','.join(str(x) for x in r) + ']' for r in rs)
body = open(sys.argv[0].replace('genallgross.py','allgross_template.m')).read()
body = body.replace('@DN@', str(DN)).replace('@DSV@', str(ds).replace(' ', '')) \
           .replace('@RS@', rows).replace('@CS@', ', '.join(cs)) \
           .replace('@TAG@', f"{D}_{N}")
open(outp, 'w').write(body)
print(f"# wrote {outp}: {len(rs)} monomials, DN={DN}")
