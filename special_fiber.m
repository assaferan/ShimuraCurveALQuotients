models := AssociativeArray();
// The points are defined as Q-schemes so it will be possible to reduce mod p 
elliptic_pts := AssociativeArray();
embeddings := AssociativeArray();
P2<x,y,z> := ProjectiveSpace(Rationals(),2);
models[<6,1,{1}>] := Curve(P2, x^2 + 3*y^2 + z^2);

models[<6,1,{1,2,3,6}>] := Curve(P2, x + 3*y + z);
for W in [{1},{1,2,3,6}] do
    elliptic_pts[<6,1,W>] := [models[<6,1,W>] meet Curve(P2,xx) : xx in [x,y,z]];
    embeddings[<6,1,W>] := map<models[<6,1,W>] -> models[<6,1,W>] | [x,y,z]>;
end for;

models[<6,1,{1,2}>] := Curve(P2, 3*x*z + y^2 + z^2);
P121<x,y,z> := ProjectiveSpace(Rationals(),[1,2,1]);
C := Curve(P121, x^2 + 3*y + z^2);
phi := map<C -> models[<6,1,{1,2}>] | [y, x*z, x^2]>;
elliptic_pts[<6,1,{1,2}>] := [phi(C meet Curve(P121,xx)) : xx in [x,y,z]];
embeddings[<6,1,{1,2}>] := phi;

_<x,y,z> := P2;
models[<6,1,{1,6}>] := Curve(P2, z^2 + 3*y^2 + x*z);
P112<x,y,z> := ProjectiveSpace(Rationals(),[1,1,2]);
C := Curve(P112, x^2 + 3*y^2 + z);
phi := map<C -> models[<6,1,{1,6}>] | [z, x*y, x^2]>;
elliptic_pts[<6,1,{1,6}>] := [phi(C meet Curve(P112,xx)) : xx in [x,y,z]];
embeddings[<6,1,{1,6}>] := phi;

_<x,y,z> := P2;
models[<6,1,{1,3}>] := Curve(P2, z^2 + 3*y^2 + x*z);
P211<x,y,z> := ProjectiveSpace(Rationals(),[2,1,1]);
C := Curve(P211, x + 3*y^2 + z^2);
phi := map<C -> models[<6,1,{1,3}>] | [x, z*y, z^2]>;
elliptic_pts[<6,1,{1,3}>] := [phi(C meet Curve(P211,xx)) : xx in [x,y,z]];
embeddings[<6,1,{1,3}>] := phi;
/*
ws := [{1,3},{1,2},{1,6}];
coeffs := [1,3,1];
for i->w in ws do
    P<x,y,z> := ProjectiveSpace(Rationals(),[1 : j in [1..i-1]] cat [2] cat [1 : j in [i+1..3]]);
    P3 := ProjectiveSpace(Rationals(),3);
    // This function might be buggy - check that it does what we want
    phi := map<P -> P3 | [m : m in MonomialsOfWeightedDegree(CoordinateRing(P),2)]>;
    C := Curve(P, &+[coeffs[j]*P.j^(2 div Gradings(P)[1][j]) : j in [1..3]]);
    models[<6,1,w>] := phi(C);
    elliptic_pts[<6,1,w>] := [phi(C meet Curve(P,xx)) : xx in [x,y,z]];
    embeddings[<6,1,w>] := phi;
end for;
*/

function get_P1_automorphism_01oo(triple)
    p0 := triple[1];
    p1 := triple[2];
    poo := triple[3];
    c0 := p1[1]*p0[2] - p0[1]*p1[2];
    d0 := poo[1]*p1[2] - p1[1]*poo[2];
    c := c0*poo[2];
    d := d0*p0[2];
    a := c0*poo[1]; // (p1-p0)*poo multiplied by all denoms
    b := p0[1]*d0; // p0(poo-p1)
    return Matrix([[a,b],[c,d]]);
end function;

function check_aut(g, origin, dest)
    g_origin := Parent(dest)![g[1,1]*origin[1] + g[1,2]*origin[2], g[2,1]*origin[1] + g[2,2]*origin[2]];
    return g_origin eq dest;
end function;

function get_P1_automorphism(P1, origin_triple, dest_triple)
    assert #origin_triple eq 3;
    g1 := get_P1_automorphism_01oo(origin_triple);
    g2 := get_P1_automorphism_01oo(dest_triple);
    ret := g2*g1^(-1);
    assert &and[check_aut(ret, origin_triple[i], dest_triple[i]) : i in [1..3]];
    return ret;
end function;

function check_aut_P1(C, C2, p, ss_pts)
    C_conic, C_to_conic := Conic(C);
    has_pt, pt := HasRationalPoint(C_conic);
    assert has_pt; // remember to handle also the case when it doesn't
    C_P1, pi := Projection(C_conic, pt);
    C2_P1 := ChangeRing(C_P1, GF(p^2));
    C2_conic := ChangeRing(C_conic, GF(p^2));
    alg_map := AlgebraMap(C_to_conic);
    polys := [alg_map(Domain(alg_map).i) : i in [1..3]];
    C2_to_conic := map<C2 -> C2_conic | polys>; 
    alg_map := AlgebraMap(pi);
    polys := [alg_map(Domain(alg_map).i) : i in [1..2]];
    conic_to_P1 := map<C2_conic -> C2_P1 | polys>;
    C2_to_P1 := C2_to_conic*conic_to_P1;
    origins := [C2_to_P1(pt) : pt in ss_pts];
    dests := [C2_to_P1(C2![Frobenius(x) : x in Eltseq(ss_pt)]) : ss_pt in ss_pts];
    if #origins lt 3 then
        num_pts := 3 - #origins;
        orig_pts := [ pt : pt in Points(C_P1) | C2_P1!Eltseq(pt) notin origins];
        origins cat:= [C2_P1!Eltseq(pt) : pt in orig_pts[1..num_pts]];
        dest_pts := [ pt : pt in Points(C_P1) | C2_P1!Eltseq(pt) notin dests];
        rat_pts := CartesianPower(dest_pts, num_pts);
        for pts in rat_pts do
            gen_diag := &or[pts[i] eq pts[j] : i,j in [1..#pts] | i lt j];
            if gen_diag then continue; end if;
            cur_dests := dests cat [C2_P1!Eltseq(pts[i]) : i in [1..num_pts]];
            g := get_P1_automorphism(C2_P1, origins, cur_dests);
            // print "Order(g) = ", Order(g);
            if Order(g) eq 2 then
                is_rational := &and[x eq Frobenius(x) : x in Eltseq(g)];
                if is_rational then return true, g; end if;
            end if;
        end for;
    else
        g := get_P1_automorphism(C2_P1, origins[1..3], dests[1..3]);
        is_rational := &and[x eq Frobenius(x) : x in Eltseq(g)];
        if is_rational then 
            if &and[check_aut(g,origins[i],dests[i]) : i in [1..#origins]] then
                return true, g;
            end if; 
        end if;
    end if;
    return false, _;
end function;


function get_cm_pts(D,N,p, W, models, elliptic_pts, embeddings)
    assert D eq 6;
    assert N eq 1;
    // This is an alternative
    // L := Universe(Eltseq(RepresentativePoint(Support(Divisor(models[<D,N>],t))[1])));
    ell_pts := [pt : pt in elliptic_pts[<D,N,W>] | IsPrime(p*Integers(NumberField(AbsolutePolynomial(L)))) where _, L := PointsOverSplittingField(pt)];
    cm_pts := [ChangeRing(pt,GF(p)) : pt in ell_pts];
    cm_pts := &cat[[*pt : pt in Points(BaseChange(cm_pt,GF(p^2)))*] : cm_pt in cm_pts];
    A := 19/24;
    B := 23/24;
    C := 3/2;
    deg := p div 24;
    _<t> := PolynomialRing(GF(p^2));
    f := &+[ &*[GF(p^2) | (A+l)*(B+l)/(C+l) : l in [0..j-1]]/Factorial(j) *t^j : j in [0..deg]];
    roots := [r[1] : r in Roots(f)];
    X := ChangeRing(models[<D,N,W>], GF(p^2));
    // the t parameter is -z/x
    // phi := ChangeRing(embeddings[<D,N,W>], GF(p^2));
    phi := embeddings[<D,N,W>];
    X0<x,y,z> := ChangeRing(Domain(phi),GF(p^2));
    alg_map := AlgebraMap(phi);
    polys := [alg_map(Domain(alg_map).i) : i in [1..Ngens(Domain(alg_map))]];
    phi_p := map<X0 -> X | polys>;
    additional_cm := &cat[[* pt : pt in Points(X meet phi_p(Scheme(X0, t0*x^Degree(z)+z^Degree(x)))) *] : t0 in roots];
    if not IsEmpty(additional_cm) then
        cm_pts cat:= additional_cm;
    end if;
    return cm_pts;
end function;


// Disprove hyperellipticity using the special fiber of X_0(D, pN) at p
function special_fiber_method(D, N, p, W, models, elliptic_pts, embeddings)
    // At the moment only works for the case of X_0(6,p) with D = 6, N = 1, p = prime
    // This is the case [FH] refers to as Figure 2, where the 
    // normalization of the special fiber X_0(D,Np)/<W, W_p>  is isomorphic to 
    // the special fiber of X_0(D,N)/W at p
    assert D eq 6;
    assert N eq 1;
    C<x,y,z> := ChangeRing(models[<D,N,W>], GF(p)); 
    C2 := BaseChange(C,GF(p^2));    
    // This assumes C is a conic
    assert IsConic(C) or IsRationalCurve(C);
    cm_pts := get_cm_pts(D,N,p, W, models, elliptic_pts, embeddings);
    ss_pts := [C2!Eltseq(pt) : pt in cm_pts | exists(x){x : x in Eltseq(pt) | x ne Frobenius(x)}];
    if IsConic(C) then
        has_pt, pt := HasRationalPoint(Conic(C));
    else
        has_pt := true;
        pt := Points(C)[1];
    end if;
    if has_pt then
        is_hyp, g := check_aut_P1(C, C2, p, ss_pts);
        return is_hyp;
    else
        auts := Automorphisms(C);
        auts2 := [map<C2 -> C2 | [AlgebraMap(a)(x),AlgebraMap(a)(y), AlgebraMap(a)(z)]> : a in auts];
        is_hyp := [&and[a2(ss_pt) eq C2![Frobenius(x) : x in Eltseq(ss_pt)] : ss_pt in ss_pts] : a2 in auts2];
        return &or is_hyp;
    end if;
end function;