# From qprime.log: does a monic cubic H in Q[X] exist with |H(0)| = N|s|, |H(2)| = N|s-2| and
# |H(-1/12) H(5/4)| = N|(s+1/12)(s-5/4)| (the norms over the three star points of discriminant d,
# read off the nine forms' values and their divisor constants)?  WITH the q-term and WITHOUT it.
# The roots of H are the Hauptmodul values at the CM points; the field it cuts out is compared with
# FieldsOfDefinitionOfCMPoint in fod.m (scratch) -- it matched at -588 and -1960 (disc -588, -1960).
from fractions import Fraction as F
import re, itertools, math, sys
log=open(sys.argv[1] if len(sys.argv)>1 else 'qprime.log').read()
C={-2:F(16,3),-1:F(320,3),9:F(20,3),10:F(100,3),11:F(1280,9),12:F(1600,9),13:F(320,3),14:F(5,4),15:F(16,3)}
def val(e):
    v=F(1)
    for sgn,coef,p in re.findall(r'([+-]?)(\d*)Log(\d+)',e):
        c=int(coef) if coef else 1
        if sgn=='-': c=-c
        v*=F(p)**c
    return v
def issq(q):
    if q<0: return False
    n,d=q.numerator,q.denominator; rn,rd=math.isqrt(n),math.isqrt(d)
    return rn*rn==n and rd*rd==d
for block in log.split('=== d = ')[1:]:
    head=block.split('\n')[0]; d=int(head.split()[0]); npts=int(re.search(r'star points (\d+)',head).group(1))
    qc=[F(x) for x in re.search(r'units of log 7\): \[ (.*?) \]',block).group(1).split(',')]
    keys=[-2,-1,9,10,11,12,13,14,15]
    vals={k:val(re.search(r'form %s\s+value (\S+)'%re.escape(str(k)),block).group(1)) for k in keys}
    for tag,shift in (("WITH",0),("WITHOUT",1)):
        v={k:vals[k]/(F(7)**(npts*qc[i]) if shift else 1) for i,k in enumerate(keys)}
        A=v[-1]/C[-1]**npts; B=v[-2]/C[-2]**npts; CE=v[14]/C[14]**npts
        print("d=%d %s: N|s|=%s N|s-2|=%s N|(s+1/12)(s-5/4)|=%s"%(d,tag,A,B,CE))
        if npts!=3: print("   (degree %d: cubic solve skipped)"%npts); continue
        found=[]
        for e0,e1,e2 in itertools.product([1,-1],repeat=3):
            c=e0*A
            def Hm(a): b=(e1*B-8-4*a-c)/2; return (F(-1,1728)+a/144-b/12+c, F(125,64)+25*a/16+5*b/4+c)
            ys=[Hm(F(t))[0]*Hm(F(t))[1]-e2*CE for t in (0,1,2)]
            r=ys[0]; p=(ys[2]-2*ys[1]+ys[0])/2; q=ys[1]-ys[0]-p
            disc=q*q-4*p*r
            if issq(disc):
                sd=F(math.isqrt(disc.numerator),math.isqrt(disc.denominator))
                for a in ((-q+sd)/(2*p),(-q-sd)/(2*p)):
                    found.append((e0,e1,e2,a,(e1*B-8-4*a-c)/2,c))
        print("   rational monic cubics H: %d"%len(found))
        for f in found: print("     signs",f[:3]," H = X^3 + (%s) X^2 + (%s) X + (%s)"%f[3:])
