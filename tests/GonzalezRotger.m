// EXTERNAL ORACLE: Gonzalez-Rotger, "Non-elliptic Shimura curves of genus one",
// J. Math. Soc. Japan 58 (2006) 927-948; arXiv:math/0612732v2.  Table 1, p.8.
//
// WHY THIS TEST EXISTS.  Until now the only external oracle in the repo was Guo-Yang, which covers
// 42 bases but NOT most of the genus-one full curves.  Eight of our committed W=[1] models had no
// published equation to check against at all -- the "10_61 evidence level": internally consistent
// and trace-formula checked, but uncorroborated.  This closes that gap for every genus-one base.
//
// The eleven (D,N) with genus(X_0(D,N)) = 1 are exactly (14,1) (15,1) (21,1) (33,1) (34,1) (46,1)
// (6,5) (6,7) (6,13) (10,3) (10,7) -- the paper's list, and independently what our own curve data
// gives.  21_1 and 33_1 have no model yet (both fail upstream), and 10_3's W=[1] key is empty.
//
// ⚠ COMPARE AN INVARIANT, NOT THE POLYNOMIAL.  Their models and ours need only be Q-EQUIVALENT
// (the paper's own relation: f2(x) = lambda^2 f1((a x+b)/(c x+d)) (c x+d)^4), so a coefficient
// diff would produce false mismatches.  We compare the Jacobian's Cremona label.
//
// ⚠ AND NOT VIA Magma's Jacobian()/EllipticCurve(): those want a rational point, and these curves
// have NONE by construction -- that is what "non-elliptic" means.  A first attempt returned ERR on
// both sides and printed a vacuous "MATCH".  Use the paper's own invariants (Section 2, p.3):
//     I = 12 a4 a0 - 3 a3 a1 + a2^2
//     J = 72 a4 a2 a0 + 9 a3 a2 a1 - 27 a4 a1^2 - 27 a3^2 a0 - 2 a2^3
//     Jac(C) : v^2 = u^3 - 27 I u - 27 J          (isomorphic over K to Jac of y^2 = f)
//
// SELF-CHECK: the paper states the Jacobian labels itself (p.8 for N=1, p.9 for N>1).  We recompute
// them from its equations and require agreement BEFORE comparing anything of ours -- so a bug in
// the I,J implementation cannot silently pass our models.
SetVerbose("ShimuraQuotients", 0);
P<x> := PolynomialRing(Rationals());

// <D, N, f(x) from Table 1, Jacobian label as STATED in the paper>
GR := [* <14,1, -x^4+13*x^2-128,               "14a2">,
         <15,1, -3*x^4-82*x^2-27,              "15a1">,
         <21,1, -7*x^4+94*x^2-343,             "21a2">,
         <33,1, -3*x^4-10*x^2-243,             "33a1">,
         <34,1, -3*x^4+26*x^3-53*x^2-26*x-3,   "34a3">,
         <46,1, -x^4+45*x^2-512,               "46a2">,
         <6,5,  -x^4+61*x^2-1024,              "30a6">,
         <6,7,  -3*x^4-34*x^2-2187,            "42a3">,
         <6,13, -x^4-115*x^2-4096,             "78a2">,
         <10,3, -2*x^4-11*x^2-32,              "30a2">,
         <10,7, -27*x^4-40*x^3+6*x^2+40*x-27,  "70a2"> *];

function jacIJ(f)
    c := [Coefficient(f,i) : i in [0..4]];
    a0:=c[1]; a1:=c[2]; a2:=c[3]; a3:=c[4]; a4:=c[5];
    I := 12*a4*a0 - 3*a3*a1 + a2^2;
    J := 72*a4*a2*a0 + 9*a3*a2*a1 - 27*a4*a1^2 - 27*a3^2*a0 - 2*a2^3;
    return MinimalModel(EllipticCurve([0,0,0,-27*I,-27*J]));
end function;

printf "Comparing genus-one models against Gonzalez-Rotger Table 1...";

// --- self-check: reproduce the paper's own stated labels from its own equations ---
for t in GR do
    got := CremonaReference(jacIJ(t[3]));
    error if got ne t[4],
        Sprintf("GR self-check FAILED at (%o,%o): recomputed Jacobian %o, paper states %o "
                * "-- the I,J implementation is wrong, so no comparison below can be trusted",
                t[1], t[2], got, t[4]);
end for;

NCMP := 0; NMISS := 0; bad := [];
for t in GR do
    D := t[1]; N := t[2]; base := Sprintf("%o_%o", D, N);
    fn := Sprintf("data/models/models_%o_%o.m", D, N);
    if not FileExists(fn) then NMISS +:= 1; continue; end if;
    models := eval (Read(fn) cat "\nreturn models;");
    key := [Integers()|1];
    if (not IsDefined(models, key)) or #models[key] eq 0 then NMISS +:= 1; continue; end if;
    e := models[key][1];
    if Type(e[2]) eq MonStgElt then NMISS +:= 1; continue; end if;   // CRV paired presentation
    fo := e[2];
    if e[3] ne 0 then fo := fo + e[3]^2/4; end if;                   // y^2+hy=f -> y^2=f+h^2/4
    if Degree(fo) gt 4 then NMISS +:= 1; continue; end if;
    NCMP +:= 1;
    ours := CremonaReference(jacIJ(fo));
    if ours ne t[4] then Append(~bad, Sprintf("%o (ours %o, GR %o)", base, ours, t[4])); end if;
end for;

error if not IsEmpty(bad),
    Sprintf("Gonzalez-Rotger oracle: %o base(s) DISAGREE with the published equation: %o",
            #bad, bad);
// ⚠ COUNT THE COMPARISONS. If models stop being found this must go red, not green-with-nothing-checked.
error if NCMP lt 8,
    Sprintf("Gonzalez-Rotger oracle: only %o comparison(s) made, expected at least 8 "
            * "(%o base(s) had no usable W=[1] model) -- something stopped being compared",
            NCMP, NMISS);
printf " ok (%o genus-one model(s) match the published equation; %o base(s) without a usable "
       * "W=[1] model)\n", NCMP, NMISS;
