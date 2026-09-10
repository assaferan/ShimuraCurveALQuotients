// tests/_crviso.m -- shared helper (leading underscore: excluded from the CI test matrix).
//
// CONSTRUCT an isomorphism between two CRV pairs, instead of calling IsIsomorphic on them.
//
// ⚠ WHY. IsIsomorphic's cost tracks PRESENTATION, not genus (tests/IsoScreen.m): a genus-7
// hyperelliptic curve settles in 0.06 s, while the genus-3 CRV pair at 14_3 runs >50 min and the
// genus-5 one at 26_3 >1 h. Four X0_*.m tests work around this by pinning a coordinate matrix
// (manual_isomorphism), which is brittle: change the pipeline's presentation -- as CMNONCOPRIME=1
// does -- and the hardcoded map stops being a map at all.
//
// ⚠ AND A SCREEN IS NOT A SUBSTITUTE. ScreenByPlaces and point counts only REFUTE: isogenous
// Jacobians agree at every prime. What follows is a PROOF: it exhibits a map and certifies it.
//
// THE IDEA. An isomorphism respecting the labelled involutions descends to the common genus-0
// base as a MOBIUS map. The y-quotient y^2 = f is hyperelliptic, where IsIsomorphic is fast, so
// take the Mobius map from there, require it to carry BOTH sides by constant squares, build the
// map, and let IsIsomorphism certify it.

// Split a CRV pair into (f, g) as one-variable polynomials, plus y's weight.
// Returns ok, f, g, wy  with  y^2 = f  and  x^2 = g  after dehomogenising the base.
function crv_split(C)
    eqs := DefiningPolynomials(C);
    if #eqs ne 2 then return false, 0, 0, 0, 0, 0, 0; end if;
    A := AmbientSpace(C);
    R := CoordinateRing(A);
    n := Rank(R);
    if n ne 4 then return false, 0, 0, 0, 0, 0, 0; end if;

    // ⚠ DO NOT ASSUME THE VARIABLE ORDER. The pipeline emits P3<x,y,s,z> with base (s,z), while
    // tests/X0_6_17.m writes P3<x,y,z,s> with base (x,s) -- hardcoding indices silently extracted
    // the wrong polynomial and HyperellipticCurve then reported "geometrically reducible".
    // Identify the roles STRUCTURALLY: a base variable occurs in BOTH equations, a fibre variable
    // (the two square roots) in exactly one.
    occ := [ [j : j in [1..n] | Degree(eqs[i], R.j) gt 0] : i in [1..2] ];
    base := [j : j in [1..n] | (j in occ[1]) and (j in occ[2])];
    fib1 := [j : j in occ[1] | j notin occ[2]];
    fib2 := [j : j in occ[2] | j notin occ[1]];
    if (#base ne 2) or (#fib1 ne 1) or (#fib2 ne 1) then return false, 0, 0, 0, 0, 0, 0; end if;

    // the y-side is the fibre variable of larger weight (its equation has the higher degree)
    wts := Gradings(A)[1];
    if wts[fib1[1]] ge wts[fib2[1]] then
        ey := eqs[1]; ex := eqs[2]; iy := fib1[1]; ix := fib2[1];
    else
        ey := eqs[2]; ex := eqs[1]; iy := fib2[1]; ix := fib1[1];
    end if;
    return true, ey, ex, base, R, iy, ix;
end function;

// Dehomogenise a CRV equation to one variable, given which base coordinate plays t.
function crv_dehom(e, base, R, which, Pt)
    t := Pt.1;
    subst := [Pt| 0 : i in [1..Rank(R)]];
    subst[base[which]] := t;
    subst[base[3-which]] := 1;
    return -Evaluate(e, subst);      // v^2 - F  ->  F
end function;

// Construct an isomorphism C -> C_ex of CRV pairs. Returns ok, psi.
function construct_crv_isomorphism(C, C_ex)
    ok1, ey1, ex1, b1, R1, iy1, ix1 := crv_split(C);
    ok2, ey2, ex2, b2, R2, iy2, ix2 := crv_split(C_ex);
    if not (ok1 and ok2) then return false, _; end if;
    Pt := PolynomialRing(Rationals());

    // ⚠ Try BOTH assignments of which base coordinate plays t on each side. The two curves need
    // not use the same convention (the pipeline and the hand-written tests do not), and a mismatch
    // would otherwise show up as a spurious "no isomorphism" rather than as the convention issue
    // it is. A wrong pairing simply fails to certify, so trying both cannot create a false pass.
    for w1 in [1,2] do
      for w2 in [1,2] do
        fo := crv_dehom(ey1, b1, R1, w1, Pt);  go := crv_dehom(ex1, b1, R1, w1, Pt);
        fg := crv_dehom(ey2, b2, R2, w2, Pt);  gg := crv_dehom(ex2, b2, R2, w2, Pt);
        if (Degree(fo) lt 3) or (Degree(fg) lt 3) then continue; end if;
        okh := true;
        try
            Cy_our := HyperellipticCurve(fo); Cy_gy := HyperellipticCurve(fg);
        catch e okh := false; end try;
        if not okh then continue; end if;
        Cy_our := HyperellipticCurve(fo); Cy_gy := HyperellipticCurve(fg);
        if Genus(Cy_our) ne Genus(Cy_gy) then continue; end if;
        if not IsIsomorphic(Cy_our, Cy_gy) then continue; end if;
        _, phi0 := IsIsomorphic(Cy_our, Cy_gy);

        cands := [phi0];
        try
            Aut, mA := AutomorphismGroup(Cy_gy);
            cands := [phi0*mA(a) : a in Aut];
        catch e ; end try;

        t := Pt.1;
        for psi in cands do
            de := DefiningEquations(psi);
            if #de lt 3 then continue; end if;
            num := Evaluate(de[1], [t, 0, 1]); den := Evaluate(de[3], [t, 0, 1]);
            if den eq 0 then continue; end if;
            mu := num/den;
            qf := Evaluate(fg, mu) / fo;  qg := Evaluate(gg, mu) / go;
            if not (IsCoercible(Rationals(), qf) and IsCoercible(Rationals(), qg)) then continue; end if;
            sf, rf := IsSquare(Rationals()!qf);  sg, rg := IsSquare(Rationals()!qg);
            if not (sf and sg) then continue; end if;

            // ⚠ THE IMAGE POLYNOMIALS LIVE IN C's RING, indexed by C_ex's coordinate POSITIONS.
            // Getting this backwards is invisible when both curves share an ambient (as in the
            // hand-check on 93_1) and wrong the moment they do not -- which is exactly the
            // pipeline-vs-testfile case, where the variable orders differ.
            SS := R1.(b1[w1]); ZZ := R1.(b1[3-w1]);
            hom_ := func<p | (p eq 0) select R1!0
                             else &+[Coefficient(p,i)*SS^i*ZZ^(Degree(p)-i) : i in [0..Degree(p)]]>;
            mn := Numerator(mu); md := Denominator(mu);
            msn := hom_(mn);
            msd := (Degree(md) eq 0) select ZZ^Degree(mn) else hom_(md);
            img := [R1| 0 : i in [1..Rank(R1)]];
            img[iy2] := rf*R1.iy1; img[ix2] := rg*R1.ix1;
            img[b2[w2]] := msn; img[b2[3-w2]] := msd;
            try
                m := map< C -> C_ex | img >;
                if IsIsomorphism(m) then return true, m; end if;
            catch e ; end try;
        end for;
      end for;
    end for;
    return false, _;
end function;
