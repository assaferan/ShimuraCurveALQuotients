// tests/_basesweep.m -- NOT a test (leading underscore: excluded from the CI matrix). Run by hand:
//     magma -b Dd:=14 Nn:=3 tests/_basesweep.m < /dev/null
//
// ⚠ WHEN A CRV PAIR WILL NOT MATCH GUO-YANG'S, SWEEP THE BASE BEFORE CONCLUDING ANYTHING ABOUT
// THE CURVE. This found base_label 5394 for 14_3 in seconds, replacing a 112-minute IsIsomorphic
// that confirmed abstract isomorphism while yielding no usable coordinate change. Same lesson as
// 26_3's 8103.
//
// Which base_label makes our W={1} CRV pair present the SAME V_4 as Guo-Yang's?
// The expensive part (Borcherds forms, CM values) runs ONCE; only the pointless-conic step
// depends on base_label, so clear the W={1} entry and replay just that step per candidate.
AttachSpec("ShimuraQuotients.spec");
import "tests/_crviso.m" : construct_crv_isomorphism;

gy := AssociativeArray();
// GY's pair, in ITS OWN weights: <weights, equations as f(a,b,c,d)>
gy[[21,2]] := <[1,3,1,1], func<a,b,c,d | [c^2 + a^2 + 3*d^2,
                  b^2 + (3*a-d)*(3*a+d)*(a^2+7*d^2)*(a^2+3*d^2)]>>;
gy[[14,3]] := <[1,2,1,1], func<a,b,c,d | [c^2 + 9*a^2 + 2*d^2, b^2 + 7*a^4 - 22*a^2*d^2 - d^4]>>;

D := StringToInteger(Dd); N := StringToInteger(Nn);
curves := GetHyperellipticCandidates();
assert exists(Xstar){X : X in curves | X`D eq D and X`N eq N and IsStarCurve(X)};
t0 := Cputime();
covers, ws := AllEquationsAboveCovers(Xstar, curves);
printf "base run: %o s\n", Cputime(t0);

assert exists(lab1){k : k in Keys(covers) | curves[k]`W eq {1}};
printf "W={1} is label %o; default bases: %o\n", lab1, Keys(covers[lab1]);

cand := {};
for oc in curves[lab1]`Covers do
    if IsDefined(covers, oc) then cand join:= Keys(covers[oc]); end if;
end for;
printf "candidate base labels: %o\n", Sort(Setseq(cand));

w, eqf := Explode(gy[[D,N]]);
Qg<a,b_,c,dd> := WeightedProjectiveSpace(Rationals(), w);
Cgy := Curve(Qg, eqf(a,b_,c,dd));
printf "GY genus %o\n", Genus(Cgy);

for bl in Sort(Setseq(cand)) do
    cp := covers;
    cp[lab1] := AssociativeArray();
    ok := true;
    try
        cp, ws2 := EquationsAbovePointlessConics(cp, ws, curves : base_label := bl);
    catch e ok := false; end try;
    if not ok or IsEmpty(Keys(cp[lab1])) then
        printf "  base %-6o : no pair produced\n", bl; continue;
    end if;
    for bb in Keys(cp[lab1]) do
        C := cp[lab1][bb];
        if Type(C) eq CrvHyp then printf "  base %-6o : hyperelliptic, not a pair\n", bl; continue; end if;
        gg := Genus(C);
        okc, psi := construct_crv_isomorphism(C, Cgy);
        printf "  base %-6o : genus %o, CONSTRUCTED ISOMORPHISM TO GY: %o\n", bl, gg, okc;
        if okc then printf "      *** eqns: %o\n", DefiningPolynomials(C); end if;
    end for;
end for;
exit;
