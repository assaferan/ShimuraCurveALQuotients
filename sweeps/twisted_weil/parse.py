import re,sys
def parse(fn):
    s=open(fn).read()
    out={}
    for blk in s.split("CreateShimuraQuot(")[1:]:
        g=lambda k: (re.search(r'<"%s", ([^>]*)>'%k,blk) or [None,None])[1]
        cid=int(g("CurveID")); out[cid]=dict(D=int(g("D")),N=int(g("N")),g=int(g("g")),sub=g("IsSubhyp"),hyp=g("IsHyp"),proof=g("TestInWhichProved"))
    return out
def parseW(fn):
    s=open(fn).read(); out={}
    for blk in s.split("CreateShimuraQuot(")[1:]:
        cid=int(re.search(r'<"CurveID", (\d+)>',blk)[1])
        w=re.search(r'<"W", \{([^}]*)\}>',blk,re.S)[1]
        w=w.replace("IntegerRing() |","")
        out[cid]=sorted(int(x) for x in re.findall(r'\d+',w))
    return out
