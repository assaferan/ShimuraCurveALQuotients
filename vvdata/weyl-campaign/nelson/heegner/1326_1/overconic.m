// Draft of the two additions to ExactDoubleCover.m, developed on a copy while a run holds the package.
//
// HyperellipticModelOverConic: for a genus-3 hyperelliptic cover whose canonical conic Q = 0 has no rational
// point, the model  y^2 = G(X1, X2, X3)  on the conic, G a quartic form: H^0(2K) = span{g_i g_j} + Q y, and y^2
// is a quartic form in the g_i modulo the conic.
//
// q_relations copied from the package (local functions are not visible from outside).
function q_relations(elts)
    coords := [[c : c in Eltseq(a), a in Eltseq(e)] : e in elts];
    den := LCM([Denominator(c) : c in cs, cs in coords]);
    polys := [[Numerator(c*den) : c in cs] : cs in coords];
    nslots := #polys[1];
    deg := Max([Degree(p) : p in ps, ps in polys] cat [0]);
    rows := [&cat[[Coefficient(p, i) : i in [0..deg]] : p in ps] : ps in polys];
    M := Matrix(Rationals(), #elts, nslots*(deg+1), rows);
    return Basis(Kernel(M));
end function;

function model_over_conic(C)
    F := C`F; gs := C`diffs;
    prods := [gs[i]*gs[j] : j in [i..3], i in [1..3]];
    idx := [[i,j] : j in [i..3], i in [1..3]];
    rel := q_relations(prods)[1];
    P2<X1,X2,X3> := ProjectiveSpace(Rationals(), 2);
    quadric := &+[rel[k]*P2.idx[k][1]*P2.idx[k][2] : k in [1..#idx]];
    // H^0(2K): functions phi with phi*omega_1^2 holomorphic
    omega := BasisOfDifferentialsFirstKind(F)[1];
    K := Divisor(omega);
    L2, m2 := RiemannRochSpace(2*K);
    assert Dimension(L2) eq 6;
    basis := [m2(L2.i) : i in [1..6]];
    // b: an element of L(2K) outside the span of the products; then the anti-invariant element
    // y = b + sum c_k q_k, determined by y^2 lying in the invariant span I_4 of the quartic monomials:
    // b^2 + 2 sum c_k (b q_k) in I_4, which is linear in the c_k.
    b := 0;
    for e in basis do
        if #q_relations(prods cat [e]) eq #q_relations(prods) then b := e; break; end if;
    end for;
    assert b ne 0;
    mons4 := [m : m in MonomialsOfDegree(CoordinateRing(P2), 4)];
    vals := [Evaluate(m, gs) : m in mons4];
    elts := [b^2] cat [b*q : q in prods] cat vals;
    rels := q_relations(elts);
    r := [v : v in rels | v[1] ne 0][1];
    cs := [r[1 + k] / (2*r[1]) : k in [1..#prods]];
    y := b + &+[cs[k]*prods[k] : k in [1..#prods]];
    rels := q_relations(vals cat [y^2]);
    r := [v : v in rels | v[#mons4 + 1] ne 0][1];
    G := -&+[r[k]*mons4[k] : k in [1..#mons4]] / r[#mons4 + 1];
    // clear denominators of G against the square of a constant in y
    den := LCM([Denominator(c) : c in Coefficients(G)]);
    G *:= den^2;  y *:= den;
    assert y^2 eq Evaluate(G, gs);
    return quadric, G, y;
end function;
