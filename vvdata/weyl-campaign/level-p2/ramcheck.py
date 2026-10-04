# Per-form check at X_0^21(2), d = -16 = 2^2 (-4): the level prime 2 divides both the conductor and d_0
# (standalone lem:ramlevel).  Input: gyforce_21_2_raw.log (DISCS, DIV, RAW lines from gyforce_21_2.m) and
# Guo-Yang's s-values (tests/_offline/GuoYang_21_2.m).  Each form is C * prod |s(d) - s(d_i)|^{m_i} over
# its divisor (the pole tau_-4 has s = oo in GY's coordinate and drops out); C is fitted at one
# fundamental point and checked at the others, then the value at -16 is compared with the code's,
# WITH the term (code) and WITHOUT (code / 2^term).  Result 2026-10-04: 9/9 with, 5/9 without
# (forms 12, 13, 14, 15 off by 2^4, 2^6, 2^8, 2^4).
import sys,re
from fractions import Fraction as F
log=open(sys.argv[1] if len(sys.argv)>1 else 'gyforce_21_2_raw.log').read().replace('\n',' ')
gyfile=sys.argv[2] if len(sys.argv)>2 else '../../../tests/_offline/GuoYang_21_2.m'
gy={}
for d,v in re.findall(r'<\s*(-\d+)\s*,\s*([^>]+)>', open(gyfile).read()):
    v=v.strip(); gy[int(d)] = None if 'Infinity' in v else F(v)
discs=[int(x) for x in re.search(r'DISCS \[(.*?)\]',log).group(1).split(',')]
divs={int(k):{int(a):int(F(b)) for a,b in re.findall(r'<\s*(-\d+),\s*([-0-9/]+)\s*>',body)} for k,body in re.findall(r'DIV (-?\d+) \[(.*?)\]',log)}
raws={int(k):[e.strip() for e in body.split(',')] for k,body in re.findall(r'RAW (-?\d+) \[(.*?)\]',log)}
term16={11:0,12:-4,13:-6,14:-8,15:-4,9:0,-1:0,10:0,-2:0}   # the term at -16 per form, from the verbose line
def val(e):
    if e in ('Logoo','Log0'): return None
    v=F(1)
    for sgn,coef,p in re.findall(r'([+-]?)([0-9/]*)Log(\d+)',e):
        c=F(coef) if coef else F(1)
        if sgn=='-': c=-c
        if c.denominator!=1: return 'irr'
        v*=F(p)**int(c)
    return v
for k in [11,12,13,14,15,9,-1,10,-2]:
    dv=divs[k]
    model=lambda d: eval('*'.join(['F(1)']+[f'abs(gy[{d}]-gy[{di}])**{mi}' for di,mi in dv.items() if gy.get(di) is not None]))
    vals={d:val(raws[k][i]) for i,d in enumerate(discs)}
    pts=[d for d in discs if gy.get(d) is not None and d not in dv and vals[d] not in (None,'irr') and d not in (-16,-100)]
    C=vals[pts[0]]/model(pts[0]); bad=[d for d in pts[1:] if vals[d]!=C*model(d)]
    truth=C*model(-16); withv=vals[-16]; without=withv/F(2)**term16[k]
    print(f"form {k:3d}: fit at {pts[0]}, checked at {len(pts)-1} points, bad={bad}; d=-16 truth={truth} with-term {'OK' if withv==truth else 'WRONG'} without-term {'OK' if without==truth else 'wrong'}")
