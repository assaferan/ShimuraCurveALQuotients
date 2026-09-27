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
// ⚠ 57_1, 82_1, 93_1 write their equations with s as the BASE coordinate and x as the CONIC
// variable, the opposite naming to 21_2/14_3 above.  In this table a is always the base, b the
// weight-(g+1) coordinate, c the conic variable, d the homogeniser -- so a = their s, c = their x.
// X_0^57(1):  y^2 = (3s+1)(3s^3+11s^2+17s+1),  x^2 = -4s^2+2s-1
gy[[57,1]] := <[1,2,1,1], func<a,b,c,d | [c^2 + 4*a^2 - 2*a*d + d^2,
                  b^2 - (3*a+d)*(3*a^3 + 11*a^2*d + 17*a*d^2 + d^3)]>>;
// X_0^82(1):  y^2 = 4s^4+4s^3+s^2-2s+1,  x^2 = -19s^2+18s-11
gy[[82,1]] := <[1,2,1,1], func<a,b,c,d | [c^2 + 19*a^2 - 18*a*d + 11*d^2,
                  b^2 - (4*a^4 + 4*a^3*d + a^2*d^2 - 2*a*d^3 + d^4)]>>;
// X_0^93(1):  y^2 = (3s^3-7s^2-3s-1)(3s^3+s^2-3s-9),  x^2 = -4s^2-6s-9
// ⚠ the JOURNAL form: arXiv v1 prints -3t for the -3s in the first factor.
gy[[93,1]] := <[1,3,1,1], func<a,b,c,d | [c^2 + 4*a^2 + 6*a*d + 9*d^2,
                  b^2 - (3*a^3 - 7*a^2*d - 3*a*d^2 - d^3)*(3*a^3 + a^2*d - 3*a*d^2 - 9*d^3)]>>;
// X_0^10(19), Guo-Yang Compositio 153 (2017) EXAMPLE 37 (body text, not the tables -- it is absent
// from A.1/A.2 because Remark 38 says this curve is not hyperelliptic over Q):
//     y^2 = -8x^6 + 57x^4 - 40x^2 + 16 ,  z^2 = 5x^2 - 32
gy[[10,19]] := <[1,3,1,1], func<a,b,c,d | [c^2 - 5*a^2 + 32*d^2,
                  b^2 + 8*a^6 - 57*a^4*d^2 + 40*a^2*d^4 - 16*d^6]>>;

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

// ⚠⚠ COMPARE THE **DEFAULT** RUN'S OWN PAIRS FIRST, and do not skip this block.  Until 2026-09-27
// this script swept only `cand` -- the bases of the OTHER cover keys -- so the pair the DEFAULT run
// actually produces was NEVER compared to Guo-Yang.  That made a negative sweep ambiguous in the
// worst way: it could not distinguish "no base presents GY's V_4" from "the default already does,
// and only the alternatives fail".  The 82_1 verdict was recorded under that ambiguity.
//   ⇒ This is also EXACTLY what the test sees.  test_AllEquationsAboveCoversSingleCurve raises its
// "matches NONE" error PER BASE ("W=%o over base %o"), so EVERY base the run produces must match --
// it is not enough that some base does.  A default base that fails is a red test.
// Free: these entries already exist from the base run above, so no replay is needed.
printf "---- pairs the DEFAULT run produced (base_label = 0) ----\n";
for b in Sort(Setseq(Keys(covers[lab1]))) do
    C := covers[lab1][b];
    if Type(C) eq CrvHyp then
        printf "  DEFAULT base %-6o : hyperelliptic, not a pair (genus %o)\n", b, Genus(C);
        continue;
    end if;
    okc := construct_crv_isomorphism(C, Cgy);
    printf "  DEFAULT base %-6o : genus %o, CONSTRUCTED ISOMORPHISM TO GY: %o\n", b, Genus(C), okc;
    printf "      weights %o\n      eqns %o\n", Gradings(Ambient(C))[1], DefiningPolynomials(C);
end for;
printf "---- alternative bases ----\n";

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
        // ⚠ Print the pair on FAILURE too.  construct_crv_isomorphism declines for two very
        // different reasons -- the pair is a different V_4, or the y-weights disagree (the 21_2
        // case) -- and without the equations a negative sweep says nothing about which, so the
        // next person pays the whole base run again to find out.
        printf "      weights %o\n      eqns %o\n",
               Gradings(Ambient(C))[1], DefiningPolynomials(C);
        printf "      GY weights %o\n", Gradings(Qg)[1];
    end for;
end for;
exit;
