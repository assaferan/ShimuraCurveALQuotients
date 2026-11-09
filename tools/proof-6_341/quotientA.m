// A = J/<t> for the 6_341 star: period lattice of J, the action of sqrt 2, the lift of t, the
// lattice of A, and whether (2 + sqrt 2) E is a principal polarisation on it.
SetColumns(0);
prec := 120;
CC := ComplexField(prec); RR := RealField(prec);
P<x> := PolynomialRing(Rationals()); Pc<z> := PolynomialRing(CC);
g := P![0,24,-47,118,37,-132,36];        // x -> 1/x of the simplified star; Weierstrass points at 0 and 3
A := AnalyticJacobian(Pc!g);
Pm := BigPeriodMatrix(A);                 // 2 x 4
printf "big period matrix computed; small period matrix symmetric error %o\n", Abs(Determinant(SmallPeriodMatrix(A) - Transpose(SmallPeriodMatrix(A))));
// endomorphisms as integer 4x4 matrices on the lattice (columns of Pm)
E := EndomorphismRing(A);
printf "endomorphism ring: rank %o\n", Dimension(E);
gens := [ ChangeRing(Matrix(b), Integers()) : b in Basis(E) ];
// find M with M^2 = 2
M := 0;
for a in [-3..3], b in [-3..3] do
    if #gens lt 2 then break; end if;
    N := a*gens[1] + b*gens[2];
    if N^2 eq 2*IdentityMatrix(Integers(), 4) then M := N; break; end if;
end for;
if Type(M) eq RngIntElt then
    // try all small combinations of all generators
    for v in CartesianPower([-2..2], #gens) do
        N := &+[ v[i]*gens[i] : i in [1..#gens] ];
        if N^2 eq 2*IdentityMatrix(Integers(), 4) then M := N; break; end if;
    end for;
end if;
require Type(M) ne RngIntElt : "no square root of 2 found among the endomorphisms";
printf "sqrt 2 on the lattice: %o\n", M;
// the symplectic form on the lattice basis: Magma's big period matrix is normalised so that the
// polarisation is the standard form J0 = [[0, I], [-I, 0]]; check Rosati symmetry of sqrt 2
J0 := Matrix(Integers(), 4, 4, [0,0,1,0, 0,0,0,1, -1,0,0,0, 0,-1,0,0]);
printf "sqrt 2 symmetric for the standard form (J0 M = M^T J0): %o\n", J0*M eq Transpose(M)*J0;
printf "   and for the other sign convention: %o\n", J0*M eq -Transpose(M)*J0;
// the lift of t = [W2 - W1]: images of the Weierstrass points x = 3 and x = 0
w1 := ToAnalyticJacobian(CC!0, CC!0, A); w2 := ToAnalyticJacobian(CC!3, CC!0, A);
v := w2 - w1;                              // 2 x 1 complex
// lattice coordinates: solve [Re Pm; Im Pm] c = [Re v; Im v]
Rm := VerticalJoin(Matrix(RR, 2, 4, [Real(e) : e in Eltseq(Pm)]), Matrix(RR, 2, 4, [Imaginary(e) : e in Eltseq(Pm)]));
rv := Matrix(RR, 4, 1, [Real(v[1,1]), Real(v[2,1]), Imaginary(v[1,1]), Imaginary(v[2,1])]);
c := Rm^-1 * rv;
printf "coordinates of t: %o\n", [RR | c[i,1] : i in [1..4]];
tau := [ Round(2*c[i,1]) : i in [1..4] ];
printf "   2*coords rounded: %o, error %o\n", tau, Maximum([Abs(2*c[i,1] - tau[i]) : i in [1..4]]);
require Maximum([Abs(2*c[i,1] - tau[i]) : i in [1..4]]) lt 10^-20 : "t is not a 2-torsion point of the lattice found";
tauv := Vector(Integers(), tau);
printf "t is 2-torsion and sqrt2 t = 0 in J: %o\n", forall{i : i in [1..4] | IsEven((tauv*Transpose(M))[i])};   // sqrt2 (tau/2) in Lambda
// Lambda' = Lambda + Z tau/2, in a basis
H := HermiteForm(Matrix(Integers(), 5, 4, &cat[[i eq j select 2 else 0 : j in [1..4]] : i in [1..4]] cat tau));   // 2 Lambda'
Bp := ChangeRing(RowSubmatrix(H, 1, 4), Rationals()) / 2;                                                           // rows: a basis of Lambda'
printf "index [Lambda' : Lambda] = %o\n", 1/Abs(Determinant(Bp));
// polarisations alpha E with alpha = a + b sqrt 2 totally positive, N(alpha) = +-2: alpha = 2 + sqrt 2 (and its conjugate 2 - sqrt 2)
for ab in [<2,1>, <2,-1>] do
    Ma := ab[1]*IdentityMatrix(Integers(), 4) + ab[2]*M;
    Gp := Bp * ChangeRing(Transpose(Ma)*J0, Rationals()) * Transpose(Bp);     // E'(x,y) = E(alpha x, y) on Lambda' (row convention)
    integral := forall{e : e in Eltseq(Gp) | IsIntegral(e)};
    printf "alpha = %o + %o sqrt2: E(alpha ., .) on Lambda' integral: %o, determinant %o, alternating %o\n", ab[1], ab[2], integral, Determinant(Gp), Gp eq -Transpose(Gp);
end for;
// the curve of (A, E') when E' is principal: symplectic basis of (Lambda', E'), small period matrix, Rosenhain, Igusa
for ab in [<2,1>, <2,-1>] do
    Ma := ab[1]*IdentityMatrix(Integers(), 4) + ab[2]*M;
    Gp := Bp * ChangeRing(Transpose(Ma)*J0, Rationals()) * Transpose(Bp);
    if not (forall{e : e in Eltseq(Gp) | IsIntegral(e)} and Abs(Determinant(Gp)) eq 1) then continue; end if;
    Gz := ChangeRing(Gp, Integers());
    Fr, Tm := FrobeniusFormAlternating(Gz);     // Tm Gz Tm^T = Fr = [[0,I],[-I,0]]
    printf "Frobenius form for alpha = %o + %o sqrt2: %o\n", ab[1], ab[2], Fr;
    Bs := ChangeRing(Tm, Rationals()) * Bp;      // rows: a symplectic basis of Lambda' for E'
    // periods of A in that basis: columns Pm * (basis vectors as column coordinate vectors)
    Pp := Pm * ChangeRing(Transpose(Bs), CC);
    O1 := Submatrix(Pp, 1, 1, 2, 2); O2 := Submatrix(Pp, 1, 3, 2, 2);
    tauA := O1^-1 * O2;
    if Imaginary(tauA[1,1]) lt 0 then tauA := O2^-1 * O1; end if;   // the other ordering if needed
    printf "small period matrix of A: symmetric to %o, Im positive definite: %o\n", Abs(tauA[1,2] - tauA[2,1]), Imaginary(tauA[1,1]) gt 0 and Imaginary(Determinant(tauA)) gt 0;
    ros := RosenhainInvariants(tauA);
    printf "Rosenhain invariants: %o\n", ros;
    f := z*(z-1)*&*[ z - r : r in ros ];
    IC := IgusaClebschInvariants(HyperellipticCurve(f));
    // normalise to absolute invariants (weights 2, 4, 6, 10)
    abs := [ IC[1]^5/IC[4], IC[2]*IC[1]^3/IC[4], IC[3]*IC[1]^2/IC[4] ];
    printf "absolute invariants: %o\n", [ Real(a) : a in abs ];
    printf "   imaginary parts: %o\n", [ Abs(Imaginary(a)) : a in abs ];
    rat := [ BestApproximation(Real(a), 10^30) : a in abs ];
    printf "   as rationals (denominator <= 10^30): %o\n", rat;
    printf "   errors: %o\n", [ Abs(Real(abs[i]) - rat[i]) : i in [1..3] ];
    // the same invariants for C itself, for comparison
    Cs := HyperellipticCurve(g); ICc := IgusaClebschInvariants(Cs);
    printf "C's absolute invariants: %o\n", [ ICc[1]^5/ICc[4], ICc[2]*ICc[1]^3/ICc[4], ICc[3]*ICc[1]^2/ICc[4] ];
    printf "RATIONAL_INVARIANTS := %o;\n", rat;
end for;
quit;
