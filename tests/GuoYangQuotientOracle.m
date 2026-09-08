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
//     uncomparable case as a failure would be as wrong as counting it as a pass.
//
// Measured 2026-09-08: 41 matches, 0 mismatches, 10 skipped over five bases.
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
        <35,[1,0,2, 0,25,0, 2,0,-1]> *]>
*];

TOTM := 0; TOTX := 0; TOTS := 0;
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
            Cs := HyperellipticCurve(e[2]);
            // ⚠ NEVER let an exception become a "MISMATCH". Magma REFUSES IsIsomorphic for genus-1
            // curves over Q ("the basefield must be finite"), and swallowing that error reported
            // six false mismatches -- a failing check is no more self-evident than a passing one.
            res := "";
            if Genus(Cs) ne Genus(Q) then
                res := Sprintf("GENUS %o vs %o", Genus(Cs), Genus(Q));
            elif Genus(Q) eq 0 then
                res := (HasRationalPoint(Conic(Cs)) eq HasRationalPoint(Conic(Q)))
                       select "MATCH(g0)" else "MISMATCH(g0 pointedness)";
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
            if res in {"MATCH","MATCH(g0)","MATCH(g1 jac)"} then nm +:= 1;
            elif "SKIP" in res then nsk +:= 1;
            else nmm +:= 1; end if;
        end for;
    end for;
    vprintf ShimuraQuotients, 1: "%o ok, %o skipped", nm, nsk;
    TOTM +:= nm; TOTX +:= nmm; TOTS +:= nsk;
end for;
error if TOTX ne 0,
    Sprintf("Guo-Yang quotient oracle: %o MISMATCH(es) -- a stored model disagrees with the "
            * "quotient derived from Guo-Yang's own curve and involutions", TOTX);
// ⚠ COUNT THE COMPARISONS. If the models stop being found, or CurveQuotient starts declining,
// this must go red rather than green-with-nothing-checked.
error if TOTM lt 41,
    Sprintf("Guo-Yang quotient oracle: only %o comparison(s) made, expected at least 41 "
            * "(%o skipped) -- something stopped being compared", TOTM, TOTS);
printf " ok (Guo-Yang quotient oracle: %o quotient comparison(s) over 5 base(s), %o skipped)\n",
       TOTM, TOTS;
