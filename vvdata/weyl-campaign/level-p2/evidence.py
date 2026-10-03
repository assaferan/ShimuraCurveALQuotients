import re
from fractions import Fraction as F
keys=[11,12,13,14,15,9,-1,10,-2]
discs=[-12,-15,-40,-60,-120,-7,-52,-28,-240,-48,-88,-132]
sgy={-7:F(1,4),-15:F(5,4),-28:F(9,4),-40:F(0),-48:F(-1,4),-52:F(1),-60:F(-1,12),-88:F(4),-120:F(2),-132:F(-1),-240:F(-25,12)}
div={-2:{-120:1},-1:{-40:1},9:{-60:1,-120:1,-15:1},10:{-60:1,-40:1,-15:1},11:{-120:1,-40:1},
     12:{-60:1,-120:1,-40:1,-15:1},13:{-40:1},14:{-60:1,-15:1},15:{-120:1}}
fib={-240:dict(zip([-2,-1,9,10,11,12,13,14,15],[2,4,2,2,4,4,4,0,2])),
     -48:dict(zip([-2,-1,9,10,11,12,13,14,15],[0,2,-4,-4,0,-4,2,-4,0]))}
def rows(fn):
    s=open(fn).read().replace('\n','')
    out={}
    for i in range(1,10):
        m=re.search(r'RAW row %d: \[(.*?)\]'%i, s); ent=[e.strip() for e in m.group(1).split(',')]
        out[keys[i-1]]=dict(zip(discs,ent))
    return out
def parse(e):
    if e in('Logoo','Log0'): return e
    if e=='0': return {}
    d={}
    for sgn,coef,p in re.findall(r'([+-]?)([0-9/]*)Log([0-9]+)',e):
        c=F(coef) if coef else F(1)
        if sgn=='-': c=-c
        d[int(p)]=c
    return d
def val(d):
    v=F(1)
    for p,c in d.items():
        assert c.denominator==1, ('irrational',d)
        v*=F(p)**int(c)
    return v
def fmt(d):
    if isinstance(d,str): return d
    return '+'.join('%sLog%d'%(c if c!=1 else '',p) for p,c in sorted(d.items())) or '0'
A=rows('pipeA.log'); B=rows('pipeB.log')
for k in [-2,-1,9,10,11,12,13,14,15]:
    # fit C at -7 from A
    a7=val(parse(A[k][-7]))
    prod=lambda d: F(1) if True else 0
    def model(d):
        v=F(1)
        for di,mi in div[k].items(): v*=abs(sgy[d]-sgy[di])**mi
        return v
    C=a7/model(-7)
    print("form key %d  C=%s" % (k,C))
    for d in [-15,-52,-28,-60,-88,-132,-120,-40,-240,-48]:
        e=parse(A[k][d])
        if isinstance(e,str) or model(d)==0: 
            print("   d=%-5d A=%s (divisor point)"%(d,A[k][d])); continue
        truth=C*model(d)
        try: av=val(e)
        except AssertionError: av='IRRATIONAL(%s)'%A[k][d]
        line="   d=%-5d truth=%-12s A=%-12s" % (d,truth,av)
        if d in fib:
            b=parse(B[k][d]); f=fib[d][k]
            bl=dict(b); bl[2]=bl.get(2,F(0))-f
            st=dict(b); st[2]=st.get(2,F(0))-2*f
            fv=val(bl); sv=val(st)
            line+=" fibre=%-12s stripped=%-12s  %s" % (fv,sv, "FIBRE OK" if fv==truth else "fibre WRONG")
            line+=("  A OK" if av==truth else "  A wrong")
        else:
            line+=("  OK" if av==truth else "  A wrong")
        print(line)
