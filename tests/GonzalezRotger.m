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
u := x;

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

// ⇒ EXHIBIT THE MAP, DO NOT MATCH AN INVARIANT.  Matching the Jacobian's Cremona label is
// NECESSARY and NOT SUFFICIENT: genus-one quartics with the same Jacobian can be INEQUIVALENT
// TORSORS of it, and the label cannot tell them apart.  That is not hypothetical here -- it is
// exactly what this upgrade found at 6_5 and 6_13 (see KNOWN_TORSOR_DRIFT below).
//
// The certificate is Gonzalez-Rotger's own relation (Section 2, p.3):
//
//     f_GR(x) = lambda^2 * (c x + d)^4 * f_ours((a x + b)/(c x + d))
//
// as an identity in Q[x], with lambda RATIONAL.  Given it, (X,Y) |-> ((aX+b)/(cX+d),
// Y/(lambda*(cX+d)^2)) is an isomorphism over Q, because Y^2 = f_GR(X) gives
// (Y/(lambda(cX+d)^2))^2 = f_GR(X)/(lambda^2 (cX+d)^4) = f_ours((aX+b)/(cX+d)).  So we check an
// exact polynomial identity and never call IsIsomorphic or Jacobian():
//   * these curves have NO rational point by construction ("non-elliptic"), so Jacobian()/
//     EllipticCurve() return ERR on both sides and print a vacuous MATCH;
//   * IsIsomorphic on a genus-0 CrvHyp is wrong on Magma 2.29-10 (Magma#125).
//
// ⚠ The lambda^2 matters.  IsGL2Equivalent decides equivalence of binary quartics MODULO ANY
// SCALAR; the curves y^2 = f are isomorphic only when that scalar is a SQUARE.  A non-square
// constant is a different torsor, which is the whole point of this check.
//
// ⚠ Asymmetry, deliberate: finding a transformation with a square constant PROVES isomorphism.
// Finding none does NOT prove non-isomorphism, since IsGL2Equivalent does not promise the full
// orbit.  So a proof is asserted, and a failure to prove is reported, never asserted upon.
function exhibit_iso(fo, fgr)
    ok, Ts := IsGL2Equivalent(fo, fgr, 4);
    if not ok then return false, _, _; end if;      // not even GL2-equivalent: no isomorphism
    R := Parent(fgr); z := R.1;
    for T in Ts do
        a, b, c, d := Explode(T);
        den := c*z + d;
        if den eq 0 then continue; end if;
        num := den^4 * Evaluate(fo, (a*z + b)/den);
        if num eq 0 or not IsCoercible(Rationals(), fgr/num) then continue; end if;
        sq, lam := IsSquare(Rationals()!(fgr/num));
        if sq then return true, T, lam; end if;
    end for;
    return false, _, _;
end function;

// ⚠ TWO COMMITTED ENTRIES ARE THE WRONG TORSOR.  At 6_5 and 6_13 the W=[1] key holds TWO entries
// that are NOT GL2-equivalent TO EACH OTHER -- genuinely different curves -- and BOTH carry the
// Jacobian label the paper states.  Only one of each pair is Gonzalez-Rotger's curve.  The
// invariant check could not see this, and since it only ever read entry [1], at 6_5 it was
// certifying the entry that is NOT the published curve.
//   6_5  : entry 1 is spurious, entry 2 is GR's curve
//   6_13 : entry 2 is spurious, entry 1 is GR's curve
// Recorded rather than asserted away: a NEW one must turn this test red.  These two entries
// should be removed from the data, which is a separate change.
KNOWN_TORSOR_DRIFT := { <6,5,1>, <6,13,2> };

NCMP := 0; NMISS := 0; NPROOF := 0; bad := []; drift := [];
for t in GR do
    D := t[1]; N := t[2]; base := Sprintf("%o_%o", D, N); fgr := t[3];
    fn := Sprintf("data/models/models_%o_%o.m", D, N);
    if not FileExists(fn) then NMISS +:= 1; continue; end if;
    models := eval (Read(fn) cat "\nreturn models;");
    key := [Integers()|1];
    if (not IsDefined(models, key)) or #models[key] eq 0 then NMISS +:= 1; continue; end if;
    proved_here := false;
    for i->e in models[key] do
        if Type(e[2]) eq MonStgElt then continue; end if;            // CRV paired presentation
        fo := e[2];
        if e[3] ne 0 then fo := fo + e[3]^2/4; end if;               // y^2+hy=f -> y^2=f+h^2/4
        if Degree(fo) gt 4 then continue; end if;
        NCMP +:= 1;
        ours := CremonaReference(jacIJ(fo));
        if ours ne t[4] then
            Append(~bad, Sprintf("%o entry %o (ours %o, GR %o)", base, i, ours, t[4]));
            continue;
        end if;
        p, T, lam := exhibit_iso(fo, fgr);
        if p then
            // certify the identity itself, not just the search's say-so
            a, b, c, d := Explode(T);
            assert fgr eq lam^2 * (c*x+d)^4 * Evaluate(fo, (a*x+b)/(c*x+d));
            NPROOF +:= 1; proved_here := true;
        elif <D,N,i> notin KNOWN_TORSOR_DRIFT then
            Append(~drift, Sprintf("%o entry %o: Jacobian %o matches but NO Q-isomorphism to the "
                                   * "published curve -- a different torsor", base, i, ours));
        end if;
    end for;
    if (not proved_here) and (#models[key] gt 0) then
        Append(~bad, Sprintf("%o: no entry could be proved isomorphic to the published curve", base));
    end if;
end for;

error if not IsEmpty(bad),
    Sprintf("Gonzalez-Rotger oracle: %o base(s)/entr(ies) DISAGREE with the published equation: %o",
            #bad, bad);
error if not IsEmpty(drift),
    Sprintf("Gonzalez-Rotger oracle: %o NEW wrong-torsor entr(ies) -- same Jacobian, not the same "
            * "curve: %o", #drift, drift);
// ⚠ COUNT THE COMPARISONS. If models stop being found this must go red, not green-with-nothing-checked.
error if NCMP lt 8,
    Sprintf("Gonzalez-Rotger oracle: only %o comparison(s) made, expected at least 8 "
            * "(%o base(s) had no usable W=[1] model) -- something stopped being compared",
            NCMP, NMISS);
// ⚠ And count the PROOFS separately: a run where every base fell back to the invariant check would
// otherwise pass silently, which is the weaker thing this section exists to stop doing.
error if NPROOF lt 11,
    Sprintf("Gonzalez-Rotger oracle: only %o exhibited isomorphism(s), expected at least 11 -- "
            * "the check has degraded to matching invariants", NPROOF);
// ---------------------------------------------------------------------------------------------
// PART 2: the AL-QUOTIENT keys.  GR give the involutions as well as the curves (table, p.8), so
// where the quartic is EVEN in x and the involution is (x,y) -> (-x,y), the quotient X/omega_m is
// simply y^2 = f(u) with u = x^2.  Comparing its Brauer class against our stored entry is what
// CAUGHT the 10_3 [1,2] drift, which three internal-consistency tests had passed for a week:
// all three committed entries were consistently WRONG, and only an oracle can arbitrate that.
//
// Class computed as ConicClasses.m does: the quaternion algebra (a, disc), a = leading coeff,
// disc = b^2-4ac.  ⚠ A LINEAR entry (a = 0) is y^2 = b u + c, which parametrises as
// u = (y^2-c)/b and is therefore RATIONAL, i.e. split -- not degenerate.  Mishandling that
// produced two false "NOT SPLIT" reports on the first run of this audit.
function conicRam(f)
    a := Coefficient(f,2); b := Coefficient(f,1); c := Coefficient(f,0);
    if a eq 0 then
        if b eq 0 then return "const"; end if;
        return [Integers()|];                      // linear: rational, split
    end if;
    return Sort(RamifiedPrimes(QuaternionAlgebra<Rationals() | a, b^2-4*a*c>));
end function;

// <D, N, m, quotient quadratic in u = x^2>  -- only bases whose GR quartic is even in x
QT := [* <14,1, 2,  -u^2+13*u-128>,    <15,1, 3,  -3*u^2-82*u-27>,
         <46,1, 2,  -u^2+45*u-512>,    <6,5,  2,  -u^2+61*u-1024>,
         <6,7,  3,  -3*u^2-34*u-2187>, <6,13, 2,  -u^2-115*u-4096>,
         <10,3, 2,  -2*u^2-11*u-32> *];
NQ := 0; qbad := [];
for t in QT do
    D:=t[1]; N:=t[2]; m:=t[3]; want := conicRam(t[4]);
    fn := Sprintf("data/models/models_%o_%o.m", D, N);
    if not FileExists(fn) then continue; end if;
    models := eval (Read(fn) cat "\nreturn models;");
    key := Sort([Integers()|1, m]);
    if (not IsDefined(models,key)) or #models[key] eq 0 then continue; end if;
    for i->e in models[key] do
        if Type(e[2]) eq MonStgElt then continue; end if;
        NQ +:= 1;
        got := conicRam(e[2]);
        if got cmpne want then
            Append(~qbad, Sprintf("%o_%o W=%o entry %o (ours ram %o, GR ram %o)",
                                  D, N, Sprint(key), i, got, want));
        end if;
    end for;
end for;
error if not IsEmpty(qbad),
    Sprintf("Gonzalez-Rotger quotient oracle: %o entry(ies) DISAGREE with the published "
            * "quotient: %o", #qbad, qbad);

// PART 3: Lemma 2.1 -- the quotient by omega_{D*N} is P^1 over Q, so W=[1,D*N] must be SPLIT.
NS := 0; sbad := [];
for t in GR do
    D:=t[1]; N:=t[2];
    fn := Sprintf("data/models/models_%o_%o.m", D, N);
    if not FileExists(fn) then continue; end if;
    models := eval (Read(fn) cat "\nreturn models;");
    key := Sort([Integers()|1, D*N]);
    if (not IsDefined(models,key)) or #models[key] eq 0 then continue; end if;
    for i->e in models[key] do
        if Type(e[2]) eq MonStgElt or e[1] ne 0 then continue; end if;
        NS +:= 1;
        r := conicRam(e[2]);
        if r cmpne [Integers()|] then
            Append(~sbad, Sprintf("%o_%o W=%o entry %o ram %o", D, N, Sprint(key), i, r));
        end if;
    end for;
end for;
error if not IsEmpty(sbad),
    Sprintf("Gonzalez-Rotger oracle: %o W=[1,D*N] entry(ies) are NOT split, contradicting "
            * "Lemma 2.1 (the quotient by omega_{D*N} is P^1 over Q): %o", #sbad, sbad);

error if NQ lt 15 or NS lt 19,
    Sprintf("Gonzalez-Rotger oracle: only %o quotient and %o splitness comparison(s) "
            * "(expected >= 15 and >= 19) -- something stopped being compared", NQ, NS);

printf " ok (%o genus-one entr(ies) checked, %o of them by an EXHIBITED isomorphism; "
       * "+ %o AL-quotient(s) + %o splitness check(s) match the published equations; "
       * "%o base(s) without a usable W=[1] model)\n", NCMP, NPROOF, NQ, NS, NMISS;
