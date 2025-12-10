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

    D0 := Curve(P2, Evaluate(f2, [X,1,Y,Z]));
    D0_con := Conic(D0);
    return D0_con, f1;
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
function maps_identify_conics(D, C)
    C0 := UnderlyingConic(C);
    D0 := UnderlyingConic(D);
    R<X,Y,Z> := CoordinateRing(Ambient(C0));
    _<x,y,z,s> := Ambient(C);
    S<a,b,c,d> := PolynomialRing(R, 4);
    zero := Evaluate(DefiningEquation(C0), [0,a*Y+b*Z,c*Y+d*Z]) - Evaluate(DefiningEquation(D0),[0,Y,Z]);
    coeffs, mons := CoefficientsAndMonomials(zero);
    Y2_eqn := &+[Coefficient(c, Y, 2)*mons[i] : i->c in coeffs];
    Z2_eqn := &+[Coefficient(c, Z, 2)*mons[i] : i->c in coeffs];
    YZ_eqn := (zero - Y2_eqn*Y^2 - Z2_eqn*Z^2) div (Y*Z);
    eqns := [Y2_eqn, Z2_eqn, YZ_eqn];
    A4<[t]> := AffineSpace(Rationals(),4);
    SS := Scheme(A4, [Evaluate(e, t) : e in eqns]);
    pts := PointSearch(SS, 10);
    return [map<D -> Ambient(C) | [x,y,pt[1]*z+pt[2]*s,pt[3]*z+pt[4]*s]> : pt in pts];
end function;