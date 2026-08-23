import sys
from fractions import Fraction as F

def parse(path):
    forms, coefs, acode, info = {}, {}, {}, None
    for line in open(path):
        t=line.split()
        if not t: continue
        if t[0]=='BASEINFO': info=(int(t[1]),int(t[2]),int(t[3]),int(t[4]))
        elif t[0]=='FORM':   forms[int(t[3])] = None if t[4]=='?' else F(t[4])
        elif t[0]=='COEF':
            k=int(t[3]); idx=(t[4], int(t[5]))
            coefs.setdefault(k,{})[idx] = coefs.setdefault(k,{}).get(idx,F(0)) + F(t[7])
        elif t[0]=='ACODE': acode[(t[3], int(t[4]))] = F(t[6])
    return info, forms, coefs, acode

def rref(M, rhs):
    M=[row[:] for row in M]; rhs=rhs[:]
    rows,cols=len(M),len(M[0]); piv=[]; r=0
    for c in range(cols):
        p=next((i for i in range(r,rows) if M[i][c]!=0), None)
        if p is None: continue
        M[r],M[p]=M[p],M[r]; rhs[r],rhs[p]=rhs[p],rhs[r]
        pv=M[r][c]; M[r]=[x/pv for x in M[r]]; rhs[r]=rhs[r]/pv
        for i in range(rows):
            if i!=r and M[i][c]!=0:
                f=M[i][c]; M[i]=[a-f*b for a,b in zip(M[i],M[r])]; rhs[i]=rhs[i]-f*rhs[r]
        piv.append(c); r+=1
        if r==rows: break
    return M,rhs,piv

path=sys.argv[1]
info, forms, coefs, acode = parse(path)
D,N,M_,det = info
idxs = sorted(acode)
print(f"X0^{D}({N})  M={M_}  det={det}")
print(f"indices ({len(idxs)}): {idxs}")
ks=[k for k in sorted(forms) if forms[k] is not None]
A=[[coefs.get(k,{}).get(i,F(0)) for i in idxs] for k in ks]
b=[forms[k] for k in ks]
R,rhs,piv = rref(A,b)
incons=[i for i in range(len(R)) if all(x==0 for x in R[i]) and rhs[i]!=0]
print(f"equations {len(ks)}, unknowns {len(idxs)}, rank {len(piv)}, free {len(idxs)-len(piv)}")
if incons:
    print("*** INCONSISTENT: the ground truth is not a linear functional of these principal parts ***")
    sys.exit()
spare = len(ks)-len(piv)
print(f"spare conditions satisfied: {spare}")
free=[c for c in range(len(idxs)) if c not in piv]
sol=[F(0)]*len(idxs)
for r,c in enumerate(piv): sol[c]=rhs[r]
print(f"\n{'index':>10} {'A_true (free vars = 0)':>24} {'A_code':>10} {'ratio':>10}")
for j,i in enumerate(idxs):
    ac=acode[i]
    ratio = 'n/a' if ac==0 else str(sol[j]/ac)
    tag = '  <-- FREE' if j in free else ''
    print(f"{str(i):>10} {str(sol[j]):>24} {str(ac):>10} {ratio:>10}{tag}")
if free:
    print(f"\nfree directions (kernel), so these A are NOT pinned by this base alone:")
    for c in free: print(f"   {idxs[c]}")
