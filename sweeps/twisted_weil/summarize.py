# Summarise results_*.out + out/*.out (deduplicated by (status, CurveID)); writes hits.txt next to this file.
import glob,collections,sys,os
T=os.path.dirname(os.path.abspath(__file__))+'/'
sys.path.insert(0,T)
from parse import parse
o=parse(T+'../../data/curves_after_D1Oracle.dat')
GL=[2555,2568,4190,5635,5639,6616,8495,7926,7932]
recs={}
for fn in sorted(glob.glob(T+'results_*.out'))+sorted(glob.glob(T+'out/*.out')):
    for l in open(fn):
        if l.startswith('RES'):
            f=l.split(); d=dict(zip(f[7::2],f[8::2])); recs[(f[1],int(f[2]))]=(f,d)
        if l.startswith('BADDIM'): print('BADDIM',l.strip())
cnt=collections.Counter(); hits=[]; ctrlfail=[]; h1fail=[]; notinv=[]; nonc=[]; m1=[]
for (st,cid),(f,d) in recs.items():
    fails=[] if d['fails']=='-' else [x.split(':') for x in d['fails'].split(';')]
    if d['notinv']!='-': notinv.append((st,cid,d['notinv']))
    if d['noncomm']!='-': nonc.append((st,cid,d['noncomm']))
    if d['minus1']!='-': m1.append((st,cid,d['minus1']))
    if st in 'HO' and fails: ctrlfail.append(f)
    if st=='U':
        if any(x[0]=='1' for x in fails): h1fail.append(cid)
        tw=[x for x in fails if x[0]!='1']
        cnt['tested', int(f[3])==1]+=1
        if tw:
            hits.append((cid,f,tw)); cnt['hit',int(f[3])==1]+=1
print('controls tested',sum(1 for k in recs if k[0] in 'HO'),'fail',len(ctrlfail))
for x in ctrlfail: print('CTRLFAIL',' '.join(x)[:300])
print('K tested',sum(1 for k in recs if k[0][0]=='K'))
print('U h=1 fails (should be none):',h1fail)
print('notinv',notinv[:10],len(notinv)); print('noncomm',nonc[:10],len(nonc)); print('minus1',m1)
print(cnt)
print('ops/p firing:',collections.Counter((x[0].split('*')[0] if '*' in x[0] else 'AL', x[1]) for _,_,tw in hits for x in tw))
print('first-fire p:',collections.Counter(min(int(x[1]) for x in tw) for _,_,tw in hits))
print('GL overlap tested:',[c for c in GL if ('U',c) in recs],'hit:',[c for c,_,_ in hits if c in GL])
d1=[c for c,f,_ in hits if int(f[3])==1]
print('D=1 hits and oracle:',[(c,o[c]['sub']) for c in d1])
with open(T+'hits.txt','w') as F:
    F.write('# CurveID D N g W  h:p:P_h(t) coeffs (all failing (h,p))  -- twisted Weil-polynomial test, P_h not in data/hypg<g>q<p>.txt\n')
    for c,f,tw in sorted(hits,key=lambda x:(int(x[1][3])*int(x[1][4]),x[0])):
        F.write('%d %s %s %s %s %s\n'%(c,f[3],f[4],f[5],f[6],';'.join(':'.join(x) for x in tw)))
print('genus of hits',collections.Counter(int(f[5]) for _,f,_ in hits), 'tested by genus',collections.Counter(int(f[5]) for (st,c),(f,d) in recs.items() if st=='U'))
