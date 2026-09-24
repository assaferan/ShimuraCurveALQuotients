// EXTERNAL ORACLE: Gonzalez-Rotger, "Non-elliptic Shimura curves of genus one",
// J. Math. Soc. Japan 58 (2006) 927-948; arXiv:math/0612732v2.
//
// ⚠ VERSION OF RECORD -- CHECKED 2026-09-23, because for Guo-Yang the journal and arXiv texts
// genuinely disagree and quoting the wrong one has already cost this project a wrong denominator
// (see [[guoyang-journal-version-differs]]). Here the risk does NOT materialise:
//   * everything below is transcribed from the arXiv v2 PDF (25 Apr 2008), which POSTDATES the
//     2006 journal version rather than preceding it -- the opposite of the Guo-Yang situation,
//     where the stale text was the arXiv one;
//   * Table 2 is IDENTICAL in v1 (Dec 2006) and v2 (Apr 2008) -- compared row by row -- so the
//     tables did not move across revisions, and the journal text sits between the two dates;
//   * Project Euclid serves the journal PDF only behind a bot check, so it was NOT read directly.
//     That residual gap is covered from the other side: every equation here is re-derived against
//     the paper's OWN stated Cremona labels before any model of ours is compared to it.
//
// THREE PUBLISHED LISTS, all transcribed here:
//   Table 1, p.8   -- 11 genus-one curves X_0(D,N)                        (PART 1, 2, 3)
//   Table 2, p.12  -- 17 genus-one AL quotients X_D^(m) = X_0(D,1)/<w_m>  (PART 4)
//   footnote 2, p.11 -- 3 further quotients that ARE elliptic over Q      (PART 5)
//
// WHY THIS TEST EXISTS.  Until now the only external oracle in the repo was Guo-Yang, which covers
// 42 bases but NOT most of the genus-one full curves.  Eight of our committed W=[1] models had no
// published equation to check against at all -- the "10_61 evidence level": internally consistent
// and trace-formula checked, but uncorroborated.  This closes that gap for every genus-one base.
//
// The eleven (D,N) with genus(X_0(D,N)) = 1 are exactly (14,1) (15,1) (21,1) (33,1) (34,1) (46,1)
// (6,5) (6,7) (6,13) (10,3) (10,7) -- the paper's list (Lemma 3.1), and independently what our own
// curve data gives.
// ⚠ This header used to add "21_1 and 33_1 have no model yet (both fail upstream), and 10_3's W=[1]
// key is empty". STALE as of 2026-09-23: all eleven now have a usable W=[1] model, and the run
// reports 0 bases without one. 21_1 was built by the HMFIT route and 33_1 by the w_1 scaling fix.
// Kept as a marked correction rather than deleted, because a comment saying a base "fails upstream"
// is exactly the kind of note that later gets quoted as a live obstruction.
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
// SCALAR; the curves y^2 = f are isomorphic only when that scalar is a SQUARE.
//
// ⚠⚠ ASYMMETRY, AND IT IS NOT A TECHNICALITY -- READ BEFORE ASSERTING ON A FAILURE.  A square
// constant PROVES isomorphism.  Finding none proves NOTHING, for two independent reasons:
//   * IsGL2Equivalent does not promise the full GL2 orbit, so the search can simply miss one; and
//   * more importantly, OUR MODEL NEED NOT BE GL2-EQUIVALENT TO THEIRS AT ALL.  A degree-2 map
//     from a genus-1 curve to P^1 is a Q-rational degree-2 divisor class, and those form a torsor
//     under E(Q).  When E(Q) is nontrivial ONE CURVE has SEVERAL INEQUIVALENT QUARTIC MODELS, and
//     E(Q) is nontrivial for every Jacobian in this table (measured: 14a2 [6], 15a1, 21a2 [2,2],
//     33a1, 34a3, 46a2, 30a6 [2,2], 42a3, 78a2 [2,2], 30a2 [2,6], 70a2 [2,2]).
// ⇒ "no exhibited isomorphism" means UNPROVED, never DISPROVED.  On 2026-09-15 this was briefly
// committed the other way round -- as "these entries are the wrong torsor" -- and retracted the
// same day.  See tests/Genus1Classes.m for the full statement of the correction.
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

// Entries for which no isomorphism to the published quartic could be exhibited.  At 6_5 and 6_13
// the W=[1] key holds two entries; one of each pair is provably GR's curve and the other is not
// provably anything.  ⚠ THIS IS NOT A DEFECT LIST.  Per the asymmetry above, the likely reading is
// that the pipeline computed the same curve over two bases which picked DIFFERENT degree-2 divisor
// classes -- both models correct, inequivalent as quartics.  Corroborating: the two entries of each
// pair have IDENTICAL everywhere-local solubility profiles (real place and every prime to 47).
//   6_5  : entry 2 is proved GR's curve; entry 1 unproved
//   6_13 : entry 1 is proved GR's curve; entry 2 unproved
// ⚠ DO NOT DELETE THESE ENTRIES on the strength of this list.  Recorded so that a NEW unproved
// entry is noticed rather than absorbed silently.
NO_EXHIBITED_ISO := { <6,5,1>, <6,13,2> };

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
        elif <D,N,i> notin NO_EXHIBITED_ISO then
            Append(~drift, Sprintf("%o entry %o: Jacobian %o matches, but no Q-isomorphism to the "
                                   * "published quartic could be EXHIBITED (unproved, not "
                                   * "disproved -- a second degree-2 class gives an inequivalent "
                                   * "model of the same curve)", base, i, ours));
        end if;
    end for;
    if (not proved_here) and (#models[key] gt 0) then
        Append(~bad, Sprintf("%o: no entry could be proved isomorphic to the published curve", base));
    end if;
end for;

error if not IsEmpty(bad),
    Sprintf("Gonzalez-Rotger oracle: %o base(s)/entr(ies) DISAGREE with the published equation: %o",
            #bad, bad);
// ⚠ This is a "something changed, look at it" guard, NOT a claim that the entry is wrong.  If the
// new entry is a legitimate second degree-2 model, add it to NO_EXHIBITED_ISO and say so.
error if not IsEmpty(drift),
    Sprintf("Gonzalez-Rotger oracle: %o entr(ies) with a matching Jacobian could not be proved "
            * "isomorphic to the published quartic: %o", #drift, drift);
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

// ---------------------------------------------------------------------------------------------
// PART 2b: THE COMPANION QUOTIENT, which is GENUS 1 and whose Jacobian GR PUBLISH.
//
// ⚠ WHY THIS EXISTS. PART 2 above checks the quotient by w_m, where m is the involution acting as
// (x,y) -> (-x,y) on an EVEN quartic: u = x^2 descends and the quotient is a CONIC, compared by its
// Brauer class. But each of those seven bases has a SECOND quotient reachable the same way, by
// w_{D*N/m} = w_{D*N} * w_m, via u = x^2 and v = x*y, giving v^2 = u*F(u). That one is GENUS ONE,
// so it carries strictly more information than a conic class -- and GR publish its Jacobian:
// p.8's table for N = 1 (the Jac(X_0(D,1)/<u.w>) column) and p.9's for N > 1 (the
// Jac(X_0(D,N)/<w_{D.N/m}>) column). Seven published values that nothing in this repo was reading.
//
// ⚠ THE COMPANION IS w_{D*N/m}, NOT some other involution: w_{D.N}*w_m = w_{D.N.m/gcd(m,D.N)^2},
// and m | D.N with gcd(m, D.N/m) = 1, so gcd(m, D.N) = m and the product is w_{D.N/m}. Checked
// against GR's own column headings, which name exactly that quotient.
//
// <D, N, companion m' = D*N/m, F(u) (the same F PART 2 uses), Jacobian AS PUBLISHED>
QT2 := [* <14,1, 7,  -u^2+13*u-128,     "14a1">,    // p.8: D=14 row, Jac(X/<u.w>) = A1
          <15,1, 5,  -3*u^2-82*u-27,    "15a2">,    // p.8: D=15 row, A2
          <46,1, 23, -u^2+45*u-512,     "46a1">,    // p.8: D=46 row, A1
          <6,5,  15, -u^2+61*u-1024,    "30a3">,    // p.9: (6,5)  m=2 -> 30A3
          <6,7,  14, -3*u^2-34*u-2187,  "42a6">,    // p.9: (6,7)  m=3 -> 42A6
          <6,13, 39, -u^2-115*u-4096,   "78a1">,    // p.9: (6,13) m=2 -> 78A1
          <10,3, 15, -2*u^2-11*u-32,    "30a1"> *]; // p.9: (10,3) m=2 -> 30A1
error if #QT2 ne #QT,
    Sprintf("GR PART 2b: %o companion row(s) against %o quotient rows -- every even-quartic base has "
            * "exactly one companion, so these must agree", #QT2, #QT);

// --- self-check: the derived companion must reproduce the PUBLISHED Jacobian, before any model of
// ours is compared to it. Same discipline as PART 1 and PART 4.
for t in QT2 do
    comp := u*(t[4]);
    got := CremonaReference(jacIJ(comp));
    error if got ne t[5],
        Sprintf("GR PART 2b self-check FAILED at %o_%o w_%o: the companion v^2 = u*F(u) has Jacobian "
                * "%o, but GR publish %o. Either F(u) is mistyped or the companion involution is "
                * "not w_{D*N/m}.", t[1], t[2], t[3], got, t[5]);
end for;

// ⚠ TWO COMPANIONS ARE A DIFFERENT TORSOR, and that is expected rather than a defect -- the same
// phenomenon as NO_EXHIBITED_ISO in PART 1. Their Jacobians agree with the publication; what
// differs is WHICH degree-2 model of that curve the pipeline stored. Recorded so a NEW one is
// noticed instead of absorbed.
NO_EXHIBITED_COMPANION := { <6,5,15>, <6,7,14> };

NC := 0; NCPROOF := 0; cbad := []; cdrift := [];
for t in QT2 do
    D := t[1]; N := t[2]; mm := t[3]; comp := P!(u*(t[4]));
    fn := Sprintf("data/models/models_%o_%o.m", D, N);
    if not FileExists(fn) then continue; end if;
    models := eval (Read(fn) cat "\nreturn models;");
    key := Sort([Integers()|1, mm]);
    if (not IsDefined(models,key)) or #models[key] eq 0 then continue; end if;
    proved := false;
    for i->e in models[key] do
        if Type(e[2]) eq MonStgElt then continue; end if;
        // ⚠ GENUS GUARD, and a negative control is what put it here. The companion is GENUS ONE by
        // construction. Pointing this row at the wrong involution lands on a genus-0 conic entry,
        // whose quartic invariants describe a SINGULAR curve, and jacIJ then dies inside Magma with
        // "Curve is singular" -- red, but pointing at EllipticCurve rather than at the mistyped row.
        // Skipping non-genus-1 entries lets the count guard below deliver that diagnosis instead.
        if e[1] ne 1 then continue; end if;
        fo := e[2];
        if e[3] ne 0 then fo := fo + e[3]^2/4; end if;
        if Degree(fo) gt 4 then continue; end if;
        NC +:= 1;
        // ⚠ CONDUCTOR BEFORE CremonaReference -- THE SAME DEFECT AS PART 4'S, AND IT RECURRED HERE
        // BECAUSE I WROTE NEW CODE WITHOUT APPLYING THE EARLIER LESSON. A wrong model throws the
        // Jacobian's conductor far outside Cremona's database range, so asking for its label first
        // dies inside Magma with "Conductor is outside the database range" -- red, but pointing at
        // the elliptic-curve database rather than at the base and key that disagree. Comparing
        // conductors first turns the commonest real failure into a message that names the entry.
        Eo := jacIJ(fo);
        Ec := jacIJ(comp);
        if Conductor(Eo) ne Conductor(Ec) then
            Append(~cbad, Sprintf("%o_%o W=[1,%o] entry %o: Jacobian has conductor %o, but the "
                                  * "published companion %o has conductor %o",
                                  D, N, mm, i, Conductor(Eo), t[5], Conductor(Ec)));
            continue;
        end if;
        ours := CremonaReference(Eo);
        if ours ne t[5] then
            Append(~cbad, Sprintf("%o_%o W=[1,%o] entry %o (ours %o, GR publish %o)",
                                  D, N, mm, i, ours, t[5]));
            continue;
        end if;
        p := exhibit_iso(fo, comp);
        if p then NCPROOF +:= 1; proved := true; end if;
    end for;
    if (not proved) and (<D,N,mm> notin NO_EXHIBITED_COMPANION) then
        Append(~cdrift, Sprintf("%o_%o W=[1,%o]: Jacobian matches the published %o, but no entry "
                                * "could be proved isomorphic to the derived companion v^2 = u*F(u) "
                                * "(unproved, NOT disproved -- a second degree-2 class gives an "
                                * "inequivalent model of the same curve)", D, N, mm, t[5]));
    end if;
end for;
error if not IsEmpty(cbad),
    Sprintf("Gonzalez-Rotger companion oracle: %o quotient(s) DISAGREE with the PUBLISHED Jacobian: "
            * "%o", #cbad, cbad);
error if not IsEmpty(cdrift),
    Sprintf("Gonzalez-Rotger companion oracle: %o quotient(s) newly unprovable against the derived "
            * "companion: %o", #cdrift, cdrift);
// ⚠ Count, as everywhere else here: seven bases all have this key, so a drop means something
// stopped being compared rather than that the data changed.
error if NC lt 7,
    Sprintf("Gonzalez-Rotger companion oracle: only %o comparison(s), expected at least 7 -- a "
            * "committed model or a W=[1,D*N/m] key stopped being found", NC);

error if NQ lt 15 or NS lt 19,
    Sprintf("Gonzalez-Rotger oracle: only %o quotient and %o splitness comparison(s) "
            * "(expected >= 15 and >= 19) -- something stopped being compared", NQ, NS);

// ---------------------------------------------------------------------------------------------
// PART 4: TABLE 2, p.12 -- the SEVENTEEN genus-one Atkin-Lehner quotients X_D^(m) = X_0(D,1)/<w_m>.
//
// ⚠ WHY THIS PART EXISTS, AND WHAT IT SAYS ABOUT THE PARTS ABOVE. Until 2026-09-23 this file
// transcribed Table 1 ONLY, so the repo's reading of "the Gonzalez-Rotger oracle" was the eleven
// genus-one FULL curves. Table 2 is a second published table, of the same kind and in the same
// paper, covering a DISJOINT set of objects -- quotient keys, not full curves -- and nothing was
// checking it. Seven of its rows had a committed model at the time it was added; all seven were
// unverified against it, on bases (39_1 55_1 62_1 69_1 77_1 94_1 178_1) whose W=[1] curve was
// already covered by Guo-Yang. ⇒ A BASE BEING "COVERED" SAID NOTHING ABOUT ITS QUOTIENT KEYS.
// This is the 69_1 lesson one level down: the sweep that wires up oracles works at base
// granularity, so a published object that is not a base is invisible to it.
//
// The quotient X_D^(m) is our key W = [1, m] at N = 1 -- the same identification PART 2 already
// makes for Table 1's quotients, and it is confirmed, not assumed: all seven comparisons match
// the published curve, and all seven do so by an EXHIBITED isomorphism (below), which a wrong
// object would not survive.
//
// ⚠ TRANSCRIPTION GUARD. Each row carries the equation AND the Cremona label of its Jacobian, read
// from the same table row, and the self-check below re-derives the label from the equation. A typo
// in either one desynchronises the pair and goes red; they would have to be wrong CONSISTENTLY to
// pass, which a slip cannot arrange. Independently, Jac(X_D^(m)) is D-new, so its conductor is
// exactly D (paper, Theorem 5.4 and Section 4) -- asserted too, which catches a row transcribed
// against the wrong D even if equation and label agree with each other.
GR2 := [* <39,13,  -7*x^4 - 24*x^3 - 34*x^2 + 24*x - 7,            "39a1">,
          <55,5,   -3*x^4 - 2*x^3 - 9*x^2 + 2*x - 3,               "55a1">,
          <62,2,   -x^4 - 8*x^3 - 78*x^2 + 248*x - 961,            "62a3">,
          <69,3,   -3*(x^4 - 156*x^2 + 6912),                        "69a2">,
          <77,11,  -11*x^4 - 19*x^2 - 16,                            "77c2">,
          <85,17,  -3*x^4 + 10*x^3 - x^2 - 10*x - 3,               "85a1">,
          <94,2,   -x^4 + 9*x^2 - 32,                                "94a2">,
          <178,89, -12*x^4 + 4*x^3 - 19*x^2 - 4*x - 12,            "178b1">,
          <210,30, -43*x^4 - 686*x^3 - 2915*x^2 - 1372*x - 172,    "210e5">,
          <210,42, -43*x^4 + 600*x^3 - 986*x^2 - 600*x - 43,       "210a6">,
          <210,70, -43*x^4 + 256*x^3 - 130*x^2 - 768*x - 387,      "210c3">,
          <210,105,-43*x^4 + 2*x^3 + 505*x^2 + 12*x - 1548,        "210b6">,
          <330,3,  -3*x^4 - 22*x^3 - 125*x^2 - 66*x - 27,          "330b3">,
          <330,22, -3*x^4 + 1358*x^2 - 177147,                       "330c3">,
          <330,33, -3*x^4 + 2846*x^2 - 2381643,                      "330d2">,
          <330,165,-3*x^4 + 614*x^2 - 19683,                         "330a2">,
          <462,154,-x^4 - 283*x^2 - 16384,                           "462b2"> *];
error if #GR2 ne 17,
    Sprintf("GR Table 2 has %o rows, not the 17 published in Lemma 4.1 -- fix the list before "
            * "trusting any verdict below", #GR2);

// --- self-check: the paper's own equations must reproduce the paper's own labels ---
// ⚠ CONDUCTOR FIRST, AND THE ORDER IS LOAD-BEARING -- it was the other way round until a negative
// control showed why that is wrong. A typo'd coefficient usually throws the conductor far out of
// Cremona's database range, so calling CremonaReference first dies inside Magma with
// "Conductor is outside the database range" -- red, but pointing at Magma's elliptic-curve
// database rather than at the mistyped row. Testing the conductor first turns the commonest real
// failure into a message that names the row and says what is wrong with it.
for t in GR2 do
    E2 := jacIJ(t[3]);
    error if Conductor(E2) ne t[1],
        Sprintf("GR Table 2 self-check FAILED at (%o,%o): Jac has conductor %o, but Jac(X_D^(m)) "
                * "is D-new and must have conductor D = %o -- either the equation is mistyped or "
                * "the row is transcribed against the wrong discriminant",
                t[1], t[2], Conductor(E2), t[1]);
    got := CremonaReference(E2);
    error if got ne t[4],
        Sprintf("GR Table 2 self-check FAILED at (%o,%o): recomputed Jacobian %o, paper states %o "
                * "-- a transcription error, or the I,J implementation is wrong", t[1], t[2], got, t[4]);
end for;

N2 := 0; N2PROOF := 0; N2MISS := 0; bad2 := []; drift2 := [];
for t in GR2 do
    D := t[1]; m := t[2]; fgr := t[3];
    fn := Sprintf("data/models/models_%o_1.m", D);
    if not FileExists(fn) then N2MISS +:= 1; continue; end if;
    models := eval (Read(fn) cat "\nreturn models;");
    key := Sort([Integers()|1, m]);
    if (not IsDefined(models, key)) or #models[key] eq 0 then N2MISS +:= 1; continue; end if;
    for i->e in models[key] do
        if Type(e[2]) eq MonStgElt then continue; end if;          // CRV paired presentation
        fo := e[2];
        if e[3] ne 0 then fo := fo + e[3]^2/4; end if;             // y^2+hy=f -> y^2=f+h^2/4
        if Degree(fo) gt 4 then continue; end if;
        N2 +:= 1;
        ours := CremonaReference(jacIJ(fo));
        if ours ne t[4] then
            Append(~bad2, Sprintf("%o_1 W=[1,%o] entry %o (ours %o, GR %o)", D, m, i, ours, t[4]));
            continue;
        end if;
        // Same asymmetry as PART 1: a square lambda PROVES it, finding none proves nothing.
        p, T, lam := exhibit_iso(fo, fgr);
        if p then
            a, b, c, d := Explode(T);
            assert fgr eq lam^2 * (c*x+d)^4 * Evaluate(fo, (a*x+b)/(c*x+d));
            N2PROOF +:= 1;
        else
            Append(~drift2, Sprintf("%o_1 W=[1,%o] entry %o: Jacobian %o matches, but no "
                                    * "Q-isomorphism to the published quartic could be EXHIBITED "
                                    * "(unproved, NOT disproved -- a second degree-2 divisor class "
                                    * "gives an inequivalent quartic model of the same curve)",
                                    D, m, i, ours));
        end if;
    end for;
end for;
error if not IsEmpty(bad2),
    Sprintf("Gonzalez-Rotger Table 2 oracle: %o quotient(s) DISAGREE with the published equation: %o",
            #bad2, bad2);
error if not IsEmpty(drift2),
    Sprintf("Gonzalez-Rotger Table 2 oracle: %o quotient(s) with a matching Jacobian could not be "
            * "proved isomorphic to the published quartic: %o", #drift2, drift2);
// ⚠ COUNT THE COMPARISONS AND THE PROOFS SEPARATELY, as PART 1 does. Ten of the seventeen rows have
// no model yet (85_1 210_1 330_1 462_1 are not built), so this cannot be raised to 17 -- but it must
// never silently fall to zero, which is what a renamed key or a moved model directory would do.
error if N2 lt 7 or N2PROOF lt 7,
    Sprintf("Gonzalez-Rotger Table 2 oracle: %o comparison(s) and %o exhibited isomorphism(s), "
            * "expected at least 7 of each (%o row(s) had no usable model) -- something stopped "
            * "being compared, or the check degraded to matching invariants", N2, N2PROOF, N2MISS);

// ---------------------------------------------------------------------------------------------
// PART 5: FOOTNOTE 2, p.11 -- three quotients that ARE elliptic curves over Q.
//
// The footnote corrects the authors' OWN earlier paper [18], which had claimed X_D^(m) fails to
// have rational points over Q_5, Q_17, Q_5 for (D,m) = (35,7), (51,3), (115,23). GR state that is
// wrong: all three have rational points everywhere locally AND a global one, coming from a CM point
// on X_D, so each is an elliptic curve over Q -- namely 35A1, 51A2 and 115A1.
//
// ⚠ A DIFFERENT CLAIM IN KIND from Table 1 and Table 2, so it needs a different check. There is no
// published quartic to compare against; the content is "this genus-one curve HAS a rational point,
// and with it is THAT elliptic curve". So we search our own model for a point and, having found
// one, require the resulting elliptic curve to be the published one. The base point does not matter
// -- different choices differ by a translation, so the isomorphism class is well defined.
FN := [* <35,7,"35a1">, <51,3,"51a2">, <115,23,"115a1"> *];
// ⚠ ONE constant, quoted by the message below -- a negative control caught this reported as a
// hardcoded "2000" while the search used a different bound, i.e. an error message that lies to the
// next reader about what was actually tried.
FN_BOUND := 2000;
NFN := 0; fnbad := [];
for t in FN do
    D := t[1]; m := t[2];
    fn := Sprintf("data/models/models_%o_1.m", D);
    if not FileExists(fn) then continue; end if;
    models := eval (Read(fn) cat "\nreturn models;");
    key := Sort([Integers()|1, m]);
    if (not IsDefined(models, key)) or #models[key] eq 0 then continue; end if;
    for i->e in models[key] do
        if Type(e[2]) eq MonStgElt then continue; end if;
        fo := e[2];
        if e[3] ne 0 then fo := fo + e[3]^2/4; end if;
        if Degree(fo) notin {3,4} then continue; end if;
        NFN +:= 1;
        // Points() needs an integral model; y^2 = f  <->  (k y)^2 = k^2 f is an isomorphism over Q
        // because the scaling is a SQUARE -- the same point PART 1's lambda^2 makes.
        k  := LCM([Denominator(c) : c in Coefficients(fo)]);
        C  := HyperellipticCurve(k^2 * fo);
        pts := Points(C : Bound := FN_BOUND);
        if #pts eq 0 then
            Append(~fnbad, Sprintf("%o_1 W=[1,%o] entry %o: GR's footnote 2 says this quotient IS "
                                   * "an elliptic curve over Q (%o), so it HAS a rational point, "
                                   * "but none was found at height bound %o. ⚠ Raise the bound "
                                   * "before concluding anything -- a search finding nothing is "
                                   * "not a proof that nothing is there", D, m, i, t[3], FN_BOUND));
            continue;
        end if;
        EE := MinimalModel(EllipticCurve(C, SetToSequence(pts)[1]));
        got := CremonaReference(EE);
        if got ne t[3] then
            Append(~fnbad, Sprintf("%o_1 W=[1,%o] entry %o: ours is %o, GR's footnote 2 states %o",
                                   D, m, i, got, t[3]));
        end if;
    end for;
end for;
error if not IsEmpty(fnbad),
    Sprintf("Gonzalez-Rotger footnote-2 oracle: %o disagreement(s): %o", #fnbad, fnbad);
// 115_1 has no model, so only two of the three can fire.
error if NFN lt 2,
    Sprintf("Gonzalez-Rotger footnote-2 oracle: only %o of the 3 published quotients were checked, "
            * "expected at least 2 (35_1 and 51_1 both have models) -- something stopped being "
            * "compared", NFN);

printf " ok (%o genus-one entr(ies) checked, %o of them by an EXHIBITED isomorphism; "
       * "+ %o AL-quotient(s) + %o companion(s) vs the PUBLISHED Jacobian (%o by exhibited iso) "
       * "+ %o splitness check(s) match the published equations; "
       * "%o base(s) without a usable W=[1] model.  Table 2: %o of 17 quotient(s) checked, all by "
       * "an exhibited isomorphism, %o row(s) not yet built.  footnote 2: %o of 3 checked)\n",
       NCMP, NPROOF, NQ, NC, NCPROOF, NS, NMISS, N2, N2MISS, NFN;
