// tests/GuoYangQuotientOracle.m -- check our stored models against EVERY Atkin-Lehner quotient
// derivable from Guo-Yang's published top curve and involutions.
//
// ⚠ WHY THIS EXISTS. tests/GuoYangEquations.m compares only the equations Guo-Yang actually PRINT,
// which for most bases is the full curve alone -- so a base with fifteen cover keys got one
// comparison. But they also print the INVOLUTIONS, and every quotient follows from those:
// CurveQuotient(AutomorphismGroup(C, [w])) is the quotient, whatever W is. That turns one
// comparison per base into one per cover key, against a source entirely external to the pipeline.
// Bases with a bespoke derivation live in their own files (GuoYangQuotients_10_19.m,
// GuoYangQuotients_22_5.m); this file is the generic sweep.
//
// ⚠ EVERY INVOLUTION IS VERIFIED to be an automorphism of the published curve before it is used,
// so a mis-transcribed table entry or a wrong composed group law fails loudly rather than
// producing a quietly wrong quotient.
//
// ⚠ THREE MAGMA TRAPS, each of which produced a WRONG VERDICT here before being handled:
//   * IsIsomorphic REFUSES genus-1 curves over Q ("the basefield must be finite"). Catching that
//     exception and reporting "mismatch" turned an unusable check into six false failures -- a
//     FAILING check is no more self-evident than a passing one. Genus 1 goes through the Jacobian
//     elliptic curve instead.
//   * Jacobian(CrvHyp of genus 1) is a JacHyp, not a curve; GenusOneModel wants a quartic;
//     EllipticCurve wants a cubic. Several routes are tried, and they compute the SAME object, so
//     this is robustness rather than shopping for a favourable answer.
//   * When no elliptic model is obtainable the case is SKIPPED, not scored. Counting an
//     uncomparable case as a failure would be as wrong as counting it as a pass. ⚠ There are now
//     ZERO skips: every one was on the GUO-YANG side, where CurveQuotient returns a plain Crv for
//     genus-1 quotients -- our own models converted fine all along. Jacobian() of a genus-1 Crv
//     yields the CrvEll directly, and all 13 recovered cases MATCH. Worth noting because a skip
//     is a silent gap: it looks like a pass in the summary line.
//
// Measured: see the count guard at the bottom for the current totals.
// X_0^15(1) contributes nothing -- its top curve has genus 1 and CurveQuotient declines there.

// GENERAL Guo-Yang oracle: derive EVERY Atkin-Lehner quotient from their published top curve and
// involutions via CurveQuotient, and compare against our stored model, label by label.
P<x> := PolynomialRing(Rationals());

// <D, N, f, [<m, image of (xx,yy,zz)> for the GENERATORS Guo-Yang print]>
// image lists are on P(1,g+1,1); e.g. (x,y)->(-1/x, y/x^4) is (xx:yy:zz) -> (-zz : yy : xx).
data := [*
  <51, 1, -(x^2+3)*(243*x^6+235*x^4-31*x^2+1),
     [* <3, [-1,0,0, 0,1,0, 0,0,1]>, <51, [1,0,0, 0,-1,0, 0,0,1]> *]>,
  <15, 1, -1/3*(x^2+3)*(x^2+243),
     [* <3, [-1,0,0, 0,1,0, 0,0,1]>, <15, [1,0,0, 0,-1,0, 0,0,1]> *]>,
  <55, 1, -(x^4-x^3+x^2+x+1)*(3*x^4+x^3-5*x^2-x+3),
     [* <5, [0,0,1, 0,1,0, -1,0,0]>, <55, [1,0,0, 0,-1,0, 0,0,1]> *]>,
  <22, 3, -27*x^8 - 308*x^6 - 2146*x^4 - 308*x^2 - 27,
     [* <2, [0,0,1, 0,-1,0, -1,0,0]>, <3, [-1,0,0, 0,1,0, 0,0,1]>,
        <66,[1,0,0, 0,-1,0, 0,0,1]> *]>,
  <15, 2, -(x^2+3)*(3*x^2+4)*(x^4-x^2+4),
     [* <2, [0,0,1, 0,-4,0, 2,0,0]>, <3, [-1,0,0, 0,1,0, 0,0,1]>,
        <5, [-1,0,0, 0,-1,0, 0,0,1]> *]>,
  <14, 5, -23*x^8 - 180*x^7 - 358*x^6 - 168*x^5 - 677*x^4 + 168*x^3 - 358*x^2 + 180*x - 23,
     [* <2, [0,0,1, 0,1,0, -1,0,0]>, <14,[1,0,0, 0,-1,0, 0,0,1]>,
        <35,[1,0,2, 0,25,0, 2,0,-1]> *]>  ,
  // ---- further level-one bases, transcribed from Table A.1 of the journal ------------------
  // Encoding: (x,y)->(-x,y) is [-1,0,0, 0,1,0, 0,0,1]; (x,y)->(x,-y) is [1,0,0, 0,-1,0, 0,0,1];
  // (x,y)->(-x,-y) is [-1,0,0, 0,-1,0, 0,0,1]; (x,y)->(-1/x, y/x^k) is [0,0,1, 0,1,0, -1,0,0];
  // (x,y)->(-1/x,-y/x^k) is [0,0,1, 0,-1,0, -1,0,0].
  <26, 1, -2*x^6 + 19*x^4 - 24*x^2 - 169,
     [* <2, [-1,0,0, 0,-1,0, 0,0,1]>, <26, [1,0,0, 0,-1,0, 0,0,1]> *]>,
  <35, 1, -(x^2+7)*(7*x^6+51*x^4+197*x^2+1),
     [* <5, [-1,0,0, 0,-1,0, 0,0,1]>, <35, [1,0,0, 0,-1,0, 0,0,1]> *]>,
  <38, 1, -16*x^6 - 59*x^4 - 82*x^2 - 19,
     [* <2, [-1,0,0, 0,-1,0, 0,0,1]>, <38, [1,0,0, 0,-1,0, 0,0,1]> *]>,
  <39, 1, -(x^4-x^3-x^2+x+1)*(7*x^4-23*x^3+5*x^2+23*x+7),
     [* <13, [0,0,1, 0,1,0, -1,0,0]>, <39, [1,0,0, 0,-1,0, 0,0,1]> *]>,
  // ⚠ Guo-Yang list w_2 and w_29 here, NOT w_58 -- transcribed as printed.
  <58, 1, -2*x^6 - 78*x^4 - 862*x^2 - 1682,
     [* <2, [-1,0,0, 0,-1,0, 0,0,1]>, <29, [1,0,0, 0,-1,0, 0,0,1]> *]>,
  <62, 1, -64*x^8 - 99*x^6 - 90*x^4 - 43*x^2 - 8,
     [* <2, [-1,0,0, 0,1,0, 0,0,1]>, <62, [1,0,0, 0,-1,0, 0,0,1]> *]>,
  <74, 1, -2*x^10 + 47*x^8 - 328*x^6 + 946*x^4 - 4158*x^2 - 1369,
     [* <2, [-1,0,0, 0,-1,0, 0,0,1]>, <74, [1,0,0, 0,-1,0, 0,0,1]> *]>,
  <86, 1, -16*x^10 + 245*x^8 - 756*x^6 - 1506*x^4 - 740*x^2 - 43,
     [* <2, [-1,0,0, 0,-1,0, 0,0,1]>, <86, [1,0,0, 0,-1,0, 0,0,1]> *]>,
  <87, 1, -(x^6-7*x^4+43*x^2+27)*(243*x^6+523*x^4+369*x^2+81),
     [* <3, [-1,0,0, 0,1,0, 0,0,1]>, <87, [1,0,0, 0,-1,0, 0,0,1]> *]>,
  <94, 1, -8*x^8 + 69*x^6 - 234*x^4 + 381*x^2 - 256,
     [* <2, [-1,0,0, 0,1,0, 0,0,1]>, <94, [1,0,0, 0,-1,0, 0,0,1]> *]>  ,
  // ---- Table A.2 (level greater than one) ---------------------------------------------------
  // ⚠ TRANSCRIPTION TRAP: the PDF text layer drops superscripts. 10_23's w_2 reads "-55 y" but is
  // -5^5 y = -3125 y (confirmed against this repo's existing ws_data for that shape), and every
  // "y/x6" is y/x^6. The automorphism check below is the safety net: a wrong coefficient makes
  // the non-trivial involutions stop preserving the curve.
  <6, 29, -64*x^12 + 813*x^10 - 3066*x^8 + 4597*x^6 - 12264*x^4 + 13008*x^2 - 4096,
     [* <2,   [-1,0,0, 0,1,0, 0,0,1]>,          // (-x, y)
        <3,   [0,0,1, 0,8,0, -2,0,0]>,          // (-2/x, 8y/x^6)  -- NOTE THE MINUS
        <174, [1,0,0, 0,-1,0, 0,0,1]> *]>,      // (x, -y)
  <6, 31, -243*x^12 + 11882*x^10 - 177701*x^8 + 803948*x^6 - 1599309*x^4 + 962442*x^2 - 177147,
     [* <2,   [0,0,1, 0,-27,0, 3,0,0]>,         // (3/x, -27y/x^6)
        <3,   [-1,0,0, 0,1,0, 0,0,1]>,          // (-x, y)
        <186, [1,0,0, 0,-1,0, 0,0,1]> *]>,      // (x, -y)
  <6, 37, -4096*x^12 - 18480*x^10 - 40200*x^8 - 51595*x^6 - 40200*x^4 - 18480*x^2 - 4096,
     [* <2,   [-1,0,0, 0,1,0, 0,0,1]>,          // (-x, y)
        <3,   [0,0,1, 0,1,0, 1,0,0]>,           // (1/x, y/x^6)
        <222, [1,0,0, 0,-1,0, 0,0,1]> *]>,      // (x, -y)
  <10, 11, -8*x^12 - 35*x^10 + 30*x^8 + 277*x^6 + 120*x^4 - 560*x^2 - 512,
     [* <10,  [0,0,1, 0,-8,0, -2,0,0]>,         // (-2/x, -8y/x^6)
        <22,  [0,0,1, 0,8,0, 2,0,0]>,           // ( 2/x,  8y/x^6)
        <110, [1,0,0, 0,-1,0, 0,0,1]> *]>,      // (x, -y)
*];

// ⚠ THIS SET IS EMPTY, AND THE STORY MATTERS. It briefly held models_87_1.m's [1,29], reported
// as a real defect on 2026-09-08. That was WRONG: the entry is <genus, f, h> for y^2 + h*y = f,
// and this test was reading only f. With h restored it matches Guo-Yang exactly. No stored model
// is known to disagree with Guo-Yang.
KNOWN_BAD := {};
// ⚠ A STORED ENTRY MAY BE <genus, f, h>, MEANING y^2 + h*y = f -- NOT y^2 = f. Nine entries
// across seven model files carry a nonzero h. Reading only e[2] silently drops it and yields a
// DIFFERENT CURVE of the same genus, which is exactly the wrong-object mistake this repo keeps
// paying for: it cost a false "defect" report against models_87_1.m on 2026-09-08, where
// 4*f + h^2 is precisely Guo-Yang's published polynomial.
function model_curve(e)
    if (#e ge 3) and (Type(e[3]) eq RngUPolElt) and (e[3] ne 0) then
        return HyperellipticCurve(e[2], e[3]);
    end if;
    return HyperellipticCurve(e[2]);
end function;

TOTM := 0; TOTX := 0; TOTS := 0; TOTKB := 0; NBASE := 0;
for d in data do
    D, N, f, gens := Explode(d);
    vprintf ShimuraQuotients, 1: "\n\tX_0^%o(%o): ", D, N;
    C := HyperellipticCurve(f);
    g := Genus(C);
    A3<xx,yy,zz> := Ambient(C);
    co := [xx,yy,zz];

    // build the maps and CHECK each is an automorphism
    mp := AssociativeArray();
    good := true;
    for gg in gens do
        m, e := Explode(gg);
        img := [ &+[e[3*(i-1)+j]*co[i] : i in [1..3]] : j in [1..3] ];
        try
            mp[m] := iso< C -> C | img, img >;
        catch err
            printf "\n  ⚠ w_%o is NOT an automorphism of Guo-Yang's published curve: %o\n",
                   m, err`Object;
            good := false;
        end try;
    end for;
    if not good then continue; end if;

    // generate the full AL group of maps by composition
    gm := [mp[k] : k in Keys(mp)];
    lab := [k : k in Keys(mp)];
    all := AssociativeArray();
    for k in Keys(mp) do all[k] := mp[k]; end for;
    for i in [1..#lab] do for j in [i+1..#lab] do
        mm := lab[i]*lab[j] div GCD(lab[i],lab[j])^2;
        if not IsDefined(all, mm) then all[mm] := gm[i]*gm[j]; end if;
    end for; end for;
    if #lab ge 3 then
        mm := lab[1]*lab[2]*lab[3] div (GCD(lab[1],lab[2])^2*GCD(lab[1]*lab[2] div GCD(lab[1],lab[2])^2, lab[3])^2);
        if not IsDefined(all, mm) then all[mm] := gm[1]*gm[2]*gm[3]; end if;
    end if;


    models := eval (Read(Sprintf("data/models/models_%o_%o.m", D, N)) cat "\nreturn models;");
    // ⚠ CHECK THE FULL CURVE TOO. For bases whose involutions are plain sign changes on an EVEN
    // polynomial, the automorphism check above cannot detect a mis-transcribed coefficient -- any
    // even f admits (x,y)->(-x,+-y). Comparing W={1} against our stored model is what actually
    // validates the transcription, and it is a new comparison in its own right for the bases
    // GuoYangEquations.m does not cover.
    ok1, es1 := IsDefined(models, [Integers()|1]);
    if ok1 and #es1 gt 0 and Type(es1[1][2]) ne MonStgElt then
        C1 := model_curve(es1[1]);
        if Genus(C1) ne g then
            printf "\n  ⚠ X_0^%o(%o) W={1}: our genus %o vs Guo-Yang's %o\n", D, N, Genus(C1), g;
            TOTX +:= 1;
        else
            o1 := false; try o1 := IsIsomorphic(C1, C); catch e ; end try;
            if o1 then TOTM +:= 1;
            else printf "\n  ⚠ X_0^%o(%o) W={1}: NOT isomorphic to Guo-Yang's published curve\n", D, N;
                 TOTX +:= 1; end if;
        end if;
    end if;
    nm := 0; nmm := 0; nsk := 0;
    for k in Sort([y : y in Keys(models)]) do
        W := {t : t in k | t ne 1};
        if IsEmpty(W) or #models[k] eq 0 then continue; end if;
        if not &and[IsDefined(all, t) : t in W] then continue; end if;
        okq := true; Q := 0;
        try
            G := AutomorphismGroup(C, [all[t] : t in W]);
            Q := CurveQuotient(G);
        catch err okq := false; end try;
        if not okq then printf "  W=%-16o CurveQuotient failed\n", Sprint(k); continue; end if;
        for e in models[k] do
            if Type(e[2]) eq MonStgElt then continue; end if;
            Cs := model_curve(e);
            // ⚠ NEVER let an exception become a "MISMATCH". Magma REFUSES IsIsomorphic for genus-1
            // curves over Q ("the basefield must be finite"), and swallowing that error reported
            // six false mismatches -- a failing check is no more self-evident than a passing one.
            res := "";
            if Genus(Cs) ne Genus(Q) then
                res := Sprintf("GENUS %o vs %o", Genus(Cs), Genus(Q));
            elif Genus(Q) eq 0 then
                // ⚠ POINTEDNESS IS NOT ISOMORPHISM. This branch used to compare only
                // HasRationalPoint on the two conics, so ANY two pointless conics "MATCHED" --
                // and 57 of the 170 comparisons (34%) land here. Over Q a conic is classified by
                // its quaternion algebra and IsIsomorphic decides it outright; negative-controlled
                // on y^2 = -x^2-1 against y^2 = -x^2-3, both POINTLESS, where the one-bit check
                // says MATCH and IsIsomorphic correctly says false.
                //
                // ⚠ BE HONEST ABOUT WHAT THIS CHANGED: nothing, today. Checked 2026-09-09, every
                // genus-0 quotient at all 20 oracle bases is a POINTED conic, and pointed conics
                // over Q are all isomorphic to P^1 -- so on the current data the old check was
                // accidentally equivalent and all 57 still MATCH. It is not equivalent in general:
                // 73 of the 281 genus-0 entries across data/models/ ARE pointless (6_5, 6_7, 6_83,
                // 82_1, 93_1), so the first oracle base with one would have degraded to a one-bit
                // check in silence. This is a guard against that, not a discovery.
                g0ok := true; C0 := 0; Q0 := 0;
                try C0 := Conic(Cs); Q0 := Conic(Q); catch e g0ok := false; end try;
                if g0ok then
                    try res := IsIsomorphic(C0, Q0) select "MATCH(g0)" else "MISMATCH(g0)";
                    catch e g0ok := false; end try;
                end if;
                // ⚠ SKIP, not a fallback to pointedness: a conic Magma will not build or compare
                // is not evidence either way, and silently degrading to the weaker check would
                // hide exactly what this change is meant to expose.
                if not g0ok then res := "SKIP(g0: no comparable conic model)"; end if;
            elif Genus(Q) eq 1 then
                // CurveQuotient returns a CrvEll here; compare JACOBIANS, both elliptic over Q
                // ⚠ Jacobian(CrvHyp of genus 1) is a JacHyp, not an elliptic curve. Go through
                // GenusOneModel on the quartic to get a CrvEll, then compare those.
                // ⚠ THE JACOBIAN ELLIPTIC CURVE IS THE INVARIANT WE WANT, but Magma reaches it
                // by different routes depending on the model: GenusOneModel wants a quartic,
                // EllipticCurve wants a cubic, and Jacobian(CrvHyp) of genus 1 returns a JacHyp
                // rather than a curve. Try each and take the first that succeeds -- they all
                // compute the same object, so this is robustness, not shopping for an answer.
                function g1ell(CC)
                    ok1 := false; E := 0;
                    if Type(CC) eq CrvEll then return true, CC; end if;
                    if Type(CC) eq CrvHyp then
                        ff := HyperellipticPolynomials(CC);
                        try E := Jacobian(GenusOneModel(ff)); ok1 := true; catch e ; end try;
                        if ok1 then return true, E; end if;
                        try E := EllipticCurve(ff); ok1 := true; catch e ; end try;
                        if ok1 then return true, E; end if;
                    end if;
                    try E := EllipticCurve(CC); ok1 := true; catch e ; end try;
                    if ok1 then return true, E; end if;
                    // ⚠ CurveQuotient returns a plain Crv for many genus-1 quotients, and that is
                    // what caused ALL the skips: our side converted fine, theirs did not.
                    // Jacobian() of a genus-1 Crv gives the CrvEll directly.
                    try E := Jacobian(CC); ok1 := Type(E) eq CrvEll; catch e ; end try;
                    if ok1 then return true, E; end if;
                    try E := Jacobian(GenusOneModel(CC)); ok1 := true; catch e ; end try;
                    if ok1 then return true, E; end if;
                    try
                        bh, CH := IsHyperelliptic(CC);
                        if bh then
                            fh := HyperellipticPolynomials(CH);
                            try E := Jacobian(GenusOneModel(fh)); ok1 := true; catch e ; end try;
                            if not ok1 then try E := EllipticCurve(fh); ok1 := true; catch e ; end try; end if;
                        end if;
                    catch e ; end try;
                    return ok1, E;
                end function;
                o1, E1 := g1ell(Cs);
                if Type(Q) eq CrvEll then o2 := true; E2 := Q; else o2, E2 := g1ell(Q); end if;
                if not (o1 and o2) then
                    res := "SKIP(g1: no elliptic model in Magma for this presentation)";
                else
                    try
                        res := IsIsomorphic(E1, E2) select "MATCH(g1 jac)" else "MISMATCH(g1 jac)";
                    catch err res := "ERROR(g1): " cat Sprint(err`Object); end try;
                end if;
            else
                try
                    res := IsIsomorphic(Cs, Q) select "MATCH" else "MISMATCH";
                catch err res := "ERROR: " cat Sprint(err`Object); end try;
            end if;
            vprintf ShimuraQuotients, 2: "\n\t  W=%o %o", Sprint(k), res;
            // ⚠ SKIP is neither a match nor a mismatch. Counting an uncomparable case as a
            // failure would be as wrong as counting it as a pass.
            isbad := <D, N, Sprint(k)> in KNOWN_BAD;
            ismatch := res in {"MATCH","MATCH(g0)","MATCH(g1 jac)"};
            if isbad then
                error if ismatch,
                    Sprintf("X_0^%o(%o) %o is listed as a KNOWN DEFECT but now MATCHES Guo-Yang -- "
                            * "the model was fixed, so remove it from KNOWN_BAD", D, N, Sprint(k));
                TOTKB +:= 1;
            elif ismatch then nm +:= 1;
            elif "SKIP" in res then nsk +:= 1;
            else nmm +:= 1; end if;
        end for;
    end for;
    vprintf ShimuraQuotients, 1: "%o ok, %o skipped", nm, nsk;
    TOTM +:= nm; TOTX +:= nmm; TOTS +:= nsk; NBASE +:= 1;
end for;
error if TOTX ne 0,
    Sprintf("Guo-Yang quotient oracle: %o MISMATCH(es) -- a stored model disagrees with the "
            * "quotient derived from Guo-Yang's own curve and involutions", TOTX);
// ⚠ COUNT THE COMPARISONS. If the models stop being found, or CurveQuotient starts declining,
// this must go red rather than green-with-nothing-checked.
error if TOTM lt 170,
    Sprintf("Guo-Yang quotient oracle: only %o comparison(s) made, expected at least 170 "
            * "(%o skipped) -- something stopped being compared", TOTM, TOTS);
printf " ok (Guo-Yang quotient oracle: %o quotient comparison(s) over %o base(s), %o skipped, "
       * "%o known defect(s) still failing)\n", TOTM, NBASE, TOTS, TOTKB;
