// get a hyperelliptic curve over a quadratic field from D

function UnderlyingConic(D)
     // assert D is given as the intersection of two hypersurfaces in P(1,2,1,1)
    P3<x,y,z,s> := Ambient(D);
    assert Gradings(P3) eq [[1,2,1,1]];
    P2<X,Y,Z> := ProjectiveSpace(Rationals(),2);
    fs := DefiningEquations(D);
    assert #fs eq 2;
    assert exists(f2){f : f in fs | Degree(f) eq 2};
    f1 := [f : f in fs | f ne f2][1];
    assert exists(y_f1_idx){i : i in [1..4] | Degree(f2, P3.i) eq 0};
    assert exists(y_f2_idx){i : i in [1..4] | Degree(f1, P3.i) eq 0};
    assert y_f1_idx ne y_f2_idx;
    P1_vars := [i : i in [1..4] | i notin [y_f1_idx, y_f2_idx]];
    assert y_f1_idx eq 2; // we need this one to be the weight 2 variable.

    // hyperelliptic variable for the conic
    assert exists(z_idx){i : i in [1..4] | Coefficient(f2, i, 2) eq 1};
    vars := [X,1,Y,Z];
    p1vars := [i : i in [1..4] | i notin [2, z_idx]];
    proj1 := [vars[z_idx],1] cat [vars[p1] : p1 in p1vars];
    proj2 := proj1[1..2] cat [proj1[4], proj1[3]];
    projs := [proj1, proj2];
    D0s := [Curve(P2, Evaluate(f2, proj)) : proj in projs];
    D0_cons := [Conic(D0) : D0 in D0s];
    return D0_cons, projs;
end function;

function HyperellipticOverQuadraticExtension(D)
    D0_con, f1 := UnderlyingConic(D);
    P2<X,Y,Z> := Ambient(D0_con);

    disc := -Discriminant(D0_con);
    K := QuadraticField(Numerator(disc)*Denominator(disc));
    D0K := ChangeRing(D0_con, K);
    P1K_to_D0K := Parametrization(D0K);
    P1K_to_D0K_alg := [AlgebraMap(P1K_to_D0K)(P2.i) : i in [1..3]];
    P2K<U,Y,V> := ProjectiveSpace(K, [1,4,1]);
    alg_xsz := [Evaluate(f, [U,V]) : f in P1K_to_D0K_alg];
    map_coords := [alg_xsz[1], Y, alg_xsz[2], alg_xsz[3]];
    HK1 := Curve(P2K,Evaluate(f1, map_coords));
    DK := ChangeRing(D, K);
    HK1_to_DK := map< HK1 -> DK | map_coords>;
    is_hyp, HK, HK1_to_HK := IsHyperelliptic(HK1);
    assert is_hyp;
    assert IsInvertible(HK1_to_DK);
    DK_to_HK := HK1_to_DK^(-1)*HK1_to_HK;
    return HK, DK_to_HK;
end function;

// This is to find maps between the conics that are the identity on the P1
// equivalent to finding an isometry between the quadratic forms
function maps_identify_conics(D, C : B := 10)
    C0s, projCs := UnderlyingConic(C);
    D0s, projDs := UnderlyingConic(D);
    C0 := C0s[1];
    projC := projCs[1];
    phis := [];
    for i->D0 in D0s do
        projD := projDs[i];
        R<X,Y,Z> := CoordinateRing(Ambient(C0));
        _<x,y,z,s> := Ambient(C);
        S<a,b,c,d> := PolynomialRing(R, 4);
        scale := Evaluate(DefiningEquation(C0),[1,0,0]) / Evaluate(DefiningEquation(D0),[1,0,0]);
        zero := Evaluate(DefiningEquation(C0), [0,a*Y+b*Z,c*Y+d*Z]) - scale * Evaluate(DefiningEquation(D0),[0,Y,Z]);
        coeffs, mons := CoefficientsAndMonomials(zero);
        Y2_eqn := &+[Coefficient(c, Y, 2)*mons[i] : i->c in coeffs];
        Z2_eqn := &+[Coefficient(c, Z, 2)*mons[i] : i->c in coeffs];
        YZ_eqn := (zero - Y2_eqn*Y^2 - Z2_eqn*Z^2) div (Y*Z);
        eqns := [Y2_eqn, Z2_eqn, YZ_eqn];
        A4<[t]> := AffineSpace(Rationals(),4);
        SS := Scheme(A4, [Evaluate(e, t) : e in eqns]);
        pts := PointSearch(SS, B);
        vars := [x,y,z,s];
        projD_inv := AssociativeArray();
        for V in [X,Y,Z] do
            projD_inv[V] := vars[Index(projD, V)];
        end for;
        valsC := AssociativeArray();
        valsC[X] := <0,0,0,0,projD_inv[X]>;
        valsC[Y] := <projD_inv[Y], projD_inv[Z], 0, 0, 0>;
        valsC[Z] := <0,0,projD_inv[Y], projD_inv[Z], 0>;
        valsC[1] := <0,0,0,0,y>;
        valsD := func<pt | [&+[pt[i]*vC[i] : i in [1..4]] + vC[5] where vC := valsC[v] : v in projC]>;
        // valsD := [(v eq X) select vars[Index(projD, X)] else (v eq Y) select pt[1]*vars[Index(projD, Y)]+pt[2]*vars[Index(projD, Y)] else pt[3]*vars[Index(projD, Z)] + pt[4]*vars[Index(projD, Z)] : v in projC | v ne 1];
        // return [map<D -> Ambient(C) | [x,y,pt[1]*z+pt[2]*s,pt[3]*z+pt[4]*s]> : pt in pts];
        phis cat:= [map<D -> Ambient(C) | valsD(pt) > : pt in pts];
    end for;
    return phis;
end function;