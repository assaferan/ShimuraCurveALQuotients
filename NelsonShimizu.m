// Evaluating the weight-2 forms of an Atkin-Lehner quotient X_0^D(N)/W at points of the upper
// half-plane by the Shimizu lift (P. Nelson, Evaluating modular forms on Shimura curves, Math.
// Comp. 84 (2015), Theorem 3.1), and what the project does with the values: Noether's count of
// quadrics through the canonical image, and for genus 3 with a marked Atkin-Lehner involution a
// hyperelliptic model.  The formula, the oldform extension, the model step, the descent from a
// quadratic field and the traps are written out in paper/nelson-method.tex.
//
// EVERYTHING HERE IS NUMERICAL: a value is a truncated lattice sum at a finite working precision,
// a count of quadrics is a reading of singular values, and a model is recognised rather than
// derived and then confirmed only by point counts, which fix the isogeny class of the Jacobian
// and not the curve.
//
// Two facts the code depends on and cannot show:
//   * Magma's QuaternionOrder(A, N) returns a different Eichler order of level N from one process
//     to the next, so a computation split over processes must rebuild ONE order from its
//     serialised Z-basis (NelsonEichlerOrderBasis, the OrderBasis parameter).
//   * The truncation error of the lattice sum falls with the cutoff only until it meets the
//     rounding error; above some ten million terms the working precision must rise with the
//     cutoff, or a residual stalls at the rounding floor and looks like a genuine small value.

declare verbose NelsonShimizu, 2;

/////////////////////////////////////////////////////////////////////////////////////////////////
// W_k by a panelled Chebyshev approximation of W_k(x) e^x.  W_k has a singularity at 0, so the
// panels are geometric: on an equal-width first panel the Chebyshev series converges only at the
// rate the Bernstein ellipse through the origin allows, which costs twenty digits.

function ChebFit(F, a, b, n, RR)
    pi := Pi(RR);
    xs := [ (a+b)/2 + (b-a)/2 * Cos((2*j-1)*pi/(2*n)) : j in [1..n] ];
    fs := [ F(x) : x in xs ];
    cs := [ RR | ];
    for m in [0..n-1] do
        s := RR ! 0;
        for j in [1..n] do s +:= fs[j] * Cos(m*(2*j-1)*pi/(2*n)); end for;
        Append(~cs, 2*s/n);
    end for;
    return cs;
end function;

function ChebEval(cs, a, b, x)
    t := (2*x - a - b)/(b - a);
    n := #cs;
    d := Parent(x) ! 0; dd := Parent(x) ! 0;
    for m in [n..2 by -1] do
        sv := d; d := 2*t*d - dd + cs[m]; dd := sv;
    end for;
    return t*d - dd + cs[1]/2;
end function;

// The derivative of a Chebyshev series, in the convention ChebEval reads: the input list stands
// for c_0/2 + sum_{m >= 1} c_m T_m(t) on [a, b], and the output is the list for the derivative
// with respect to x, so ChebEval applied to it returns d/dx of the original.  Differentiating the
// fit is exact for the fit, which is what is wanted: the kernel enters the lattice sum only
// through this approximation, so its derivative must be the derivative of the same function.
function ChebDeriv(cs, a, b)
    n := #cs;
    RR := Parent(cs[1]);
    A := [ RR | cs[1]/2 ] cat [ cs[j+1] : j in [1..n-1] ];
    Bc := [ RR | 0 : j in [1..n+1] ];
    for j in [n-1 .. 1 by -1] do
        Bc[j] := Bc[j+2] + 2*j*A[j+1];
    end for;
    // The recurrence already returns the constant term in the halved convention ChebEval reads --
    // for f = T_1 it gives B_0 = 2 where f' = 1 -- so it is passed through unscaled.  Doubling it
    // here is the one error this function invites, and it is wrong by a bounded amount rather than
    // grossly, so it survives a casual look.
    sc := 2/(b - a);
    return [ RR | Bc[1]*sc ] cat [ Bc[j+1]*sc : j in [1..n-1] ];
end function;

function WkDirect(k, x, RR)
    eps := RR ! 10^(-(Precision(RR)+5));
    s := RR ! 0; n := 1;
    while n le 400 do
        u := n*x;
        t := (RR ! n)^k * (u*KBessel(k-1, u) - KBessel(k, u));
        s +:= t;
        if n ge 2 and Abs(t) lt eps*(1 + Abs(s)) then break; end if;
        n +:= 1;
    end while;
    return s / 2^(k-1);
end function;

function WkTable(k, cutmax, RR)
    a := RR ! 0.4; b := RR ! (cutmax + 1); r := RR ! 1.5;
    nps := Maximum(1, Ceiling(Log(b/a)/Log(r)));
    edges := [ a * r^i : i in [0..nps] ];
    F := func< x | WkDirect(k, x, RR) * Exp(x) >;
    coeffs := [ ChebFit(F, edges[i], edges[i+1], Precision(RR) + 14, RR) : i in [1..nps] ];
    dcoeffs := [ ChebDeriv(coeffs[i], edges[i], edges[i+1]) : i in [1..nps] ];
    return <edges, coeffs, Log(r), a, nps, dcoeffs>;
end function;

function WkFast(P, x)
    edges := P[1]; coeffs := P[2]; lr := P[3]; a := P[4]; nps := P[5];
    i := Minimum(nps, Maximum(1, Ceiling(Log(x/a)/lr)));
    return ChebEval(coeffs[i], edges[i], edges[i+1], x) * Exp(-x);
end function;

// The kernel and its derivative together.  The table fits F(x) = W_k(x) e^x, so W_k = C e^{-x}
// and W_k' = (C' - C) e^{-x}; the two Chebyshev evaluations share the panel lookup and the
// exponential, which is why asking for the derivative costs little beyond the second series.
function WkFastBoth(P, x)
    edges := P[1]; coeffs := P[2]; lr := P[3]; a := P[4]; nps := P[5]; dcoeffs := P[6];
    i := Minimum(nps, Maximum(1, Ceiling(Log(x/a)/lr)));
    e := Exp(-x);
    C := ChebEval(coeffs[i], edges[i], edges[i+1], x);
    Cd := ChebEval(dcoeffs[i], edges[i], edges[i+1], x);
    return C*e, (Cd - C)*e;
end function;

/////////////////////////////////////////////////////////////////////////////////////////////////
// Orders

intrinsic NelsonEichlerOrderBasis(D::RngIntElt, N::RngIntElt) -> SeqEnum
{A Z-basis of an Eichler order of level N in QuaternionAlgebra(D), as coordinate vectors in the
 standard basis 1, i, j, k of that algebra.  QuaternionOrder(A, N) is not deterministic across
 processes, so a computation that is split over processes, or that must be reproduced later,
 fixes its order through this sequence and NelsonEichlerOrder.}
    A := QuaternionAlgebra(D); O := QuaternionOrder(A, N);
    require Discriminant(O) eq D*N : "the order does not have discriminant D*N";
    return [ Eltseq(A ! x) : x in ZBasis(O) ];
end intrinsic;

intrinsic NelsonEichlerOrder(D::RngIntElt, N::RngIntElt, basis::SeqEnum) -> AlgQuatOrd
{The Eichler order of level N in QuaternionAlgebra(D) with the given Z-basis (coordinate vectors
 in the standard basis), as written by NelsonEichlerOrderBasis.}
    A := QuaternionAlgebra(D);
    O := QuaternionOrder([ A ! v : v in basis ]);
    require Discriminant(O) eq D*N : "the rebuilt order does not have discriminant D*N";
    require forall{ v : v in basis | IsCoercible(O, A ! v) } : "the basis does not lie in the rebuilt order";
    return O;
end intrinsic;

// the Eichler order of level M between O_N and the maximal order O_1 containing it
function SuperOrder(ON, O1, M)
    A := Algebra(ON);
    return QuaternionOrder([ A ! x : x in Basis(ON) ] cat [ A ! (M*x) : x in Basis(O1) ]);
end function;

// a basis of J_m = {x in O_N : trd(x y) in m Z for all y in O_N}, the two-sided ideal of norm m
function TraceIdealBasis(ON, m)
    A := Algebra(ON);
    bas := [ A ! x : x in Basis(ON) ];
    T := Matrix(Integers(), 4, 4, [ Trace(bas[i]*bas[j]) : i, j in [1..4] ]);
    K := Kernel(Matrix(Integers(m), T));
    rows := [ [ Integers() ! c : c in Eltseq(K.i) ] : i in [1..Dimension(K)] ]
            cat [ Eltseq(m*IdentityMatrix(Integers(), 4)[i]) : i in [1..4] ];
    B := BasisMatrix(Lattice(Matrix(Integers(), rows)));
    assert Abs(Determinant(B)) eq m^2;
    return [ &+[ B[i][j]*bas[j] : j in [1..4] ] : i in [1..4] ];
end function;

// an element of O_N of reduced norm exactly m normalising O_N: a short vector of J_m for the
// Frobenius majorant of the real splitting.  Magma's Normalizer is for definite orders only.
function ALElement(ON, m, mrep)
    A := Algebra(ON);
    jb := TraceIdealBasis(ON, m);
    R30 := RealField(30);
    Ms := [ ChangeRing(mrep(b), R30) : b in jb ];
    G := Matrix(Integers(), 4, 4,
        [ Round(10^6 * (&+[ Ms[a][r][c]*Ms[b][r][c] : r, c in [1..2] ])/2) : a, b in [1..4] ]);
    for a in [1..4] do for b in [a+1..4] do G[b][a] := G[a][b]; end for; G[a][a] +:= 1; end for;
    L := LatticeWithGram(G);
    bound := 10^6 * m * 4;
    while bound le 10^6 * m * 10^8 do
        for v in ShortVectors(L, bound) do
            x := &+[ (Integers() ! v[1][i]) * jb[i] : i in [1..4] ];
            if Norm(x) ne m then continue; end if;
            if forall{ b : b in Basis(ON) | x*(A!b)*x^-1 in ON } then return x; end if;
        end for;
        bound *:= 4;
    end while;
    error "no Atkin-Lehner element of norm", m;
end function;

function ALMul(u, v) gg := GCD(u, v); return (u div gg) * (v div gg); end function;

function ALGenerators(W)
    span := {Integers() | 1}; gens := [];
    for m in Sort(SetToSequence(W)) do
        if m in span then continue; end if;
        Append(~gens, m);
        span := { ALMul(s, m) : s in span } join span;
    end for;
    return gens;
end function;

/////////////////////////////////////////////////////////////////////////////////////////////////
// The evaluator at one level

NelLevel := recformat< D, N, k, A, O, mrep, DN, ds, latb, latQ, prec, Wtab, cutmax, scl >;

function TraceGram(bas, A)
    n := #bas;
    return Matrix(Rationals(), n, n, [ Trace(A ! (bas[i]*bas[j])) : i, j in [1..n] ]);
end function;

function NormGram(bas, A)
    n := #bas;
    return Matrix(Rationals(), n, n,
        [ Trace(A ! (bas[i]*Conjugate(A ! bas[j])))/2 : i, j in [1..n] ]);
end function;

// R^(d) = R + (DN/d) R^dual as a lattice, not as a Hermite form of the stacked generators, which
// can return a basis of the right determinant spanning the wrong lattice
function RdBasis(T, DN, d)
    n := Nrows(T);
    gens := VerticalJoin(IdentityMatrix(Rationals(), n), (DN div d)*T^(-1));
    den := LCM([ Denominator(x) : x in Eltseq(gens) ]);
    return ChangeRing(BasisMatrix(Lattice(ChangeRing(den*gens, Integers()))), Rationals())/den;
end function;

function SigZ(z, CC)
    x := Real(z); y := Imaginary(z); r := Sqrt(y);
    return Matrix(CC, 2, 2, [ r, x/r, 0, 1/r ]);
end function;

function SigZi(z, CC)
    x := Real(z); y := Imaginary(z); r := Sqrt(y);
    return Matrix(CC, 2, 2, [ 1/r, -x/r, 0, r ]);
end function;

// d/dz of SigZ.  With x = (z + conj z)/2 and y = (z - conj z)/2i one has d/dz y = -i/2, hence
// d/dz sqrt(y) = -i/(4 sqrt y), d/dz (1/sqrt y) = i/(4 y^{3/2}) and
// d/dz (x/sqrt y) = 1/(2 sqrt y) + i x/(4 y^{3/2}).  This is the only place z enters the lattice
// sum, so it is the only derivative the chain rule needs.
function SigZd(z, CC)
    x := Real(z); y := Imaginary(z); r := Sqrt(y); i_ := CC.1;
    return Matrix(CC, 2, 2, [ -i_/(4*r), 1/(2*r) + i_*x/(4*y*r), 0, i_/(4*y*r) ]);
end function;

function LevelSetup(D, N, k, prec, cutmax, Wtab, O, mrep)
    A := Algebra(O);
    assert Discriminant(O) eq D*N;
    bas := [ A ! b : b in ZBasis(O) ];
    T := TraceGram(bas, A);
    DN := D*N;
    latb := AssociativeArray(); latQ := AssociativeArray();
    for d in Divisors(DN) do
        B := RdBasis(T, DN, d);
        bd := [ &+[ B[i][j]*bas[j] : j in [1..4] ] : i in [1..4] ];
        Q := NormGram(bd, A);
        // the reduced norms on R^(d) have denominators dividing d; a wrong lattice shows here
        assert forall{ i : i in [1..4] | IsIntegral(d*Q[i][i]) };
        latb[d] := bd; latQ[d] := Q;
    end for;
    return rec< NelLevel | D := D, N := N, k := k, A := A, O := O, mrep := mrep, DN := DN,
                           ds := Divisors(DN), latb := latb, latQ := latQ, prec := prec,
                           Wtab := Wtab, cutmax := cutmax, scl := 10^16 >;
end function;

// Theorem 3.1 for several eigenforms in one enumeration: afs[i][n] = a_f(n) of the i-th form.
// Deriv also returns d/dz2 of each value, as a third return value.  The point of doing it here
// rather than by a difference quotient outside is that the enumeration of short vectors -- the
// whole cost -- is shared: asking for the derivative adds a second Chebyshev evaluation and a few
// products per term, against the dozen extra evaluations a Cauchy circle would need.
//
// Each term is not holomorphic in z2, only the sum is; but d/dz may still be applied term by term,
// since it is the d/dz-bar parts that cancel and not the d/dz ones.
function EvalMulti(dat, z1, z2, afs, cut, CC : Deriv := false)
    k := dat`k; RR := RealField(Precision(CC)); pi := Pi(RR); i_ := CC.1;
    assert IsEven(k) and cut le dat`cutmax;
    nf := #afs; nmax := #afs[1];
    s1i := SigZi(z1, CC); s2 := SigZ(z2, CC);
    s2d := Deriv select SigZd(z2, CC) else s2;
    C := cut/(4*pi); Pb := 2*C^2;
    tot := [ CC | 0 : i in [1..nf] ];
    dtot := [ CC | 0 : i in [1..nf] ];
    nterms := 0;
    // the Gram matrix is rounded to decide which vectors the enumeration offers; each one is then
    // accepted or rejected by 4 pi |X| <= cut at full precision, so the rounding costs no accuracy
    SCL := dat`scl;
    for d in dat`ds do
        bs := dat`latb[d]; Qn := dat`latQ[d];
        Ms := [ s1i * ChangeRing(dat`mrep(b), CC) * s2 : b in bs ];
        Xs := [ ((M[1][1]+M[2][2]) + i_*(M[1][2]-M[2][1]))/2 : M in Ms ];
        // M has REAL entries, so X and its conjugate are the two functionals below applied to M;
        // applying the same two to dM/dz gives dX/dz and d(conj X)/dz.  The Gram matrix, and hence
        // the enumeration, is built from M alone and is untouched by this.
        Xzs := [ CC | ]; Xbzs := [ CC | ];
        if Deriv then
            Ps := [ s1i * ChangeRing(dat`mrep(b), CC) * s2d : b in bs ];
            Xzs  := [ ((P[1][1]+P[2][2]) + i_*(P[1][2]-P[2][1]))/2 : P in Ps ];
            Xbzs := [ ((P[1][1]+P[2][2]) - i_*(P[1][2]-P[2][1]))/2 : P in Ps ];
        end if;
        G := Matrix(RR, 4, 4, [ RR ! 0 : j in [1..16] ]);
        for a in [1..4] do for b in [a..4] do
            v := (&+[ Real(Ms[a][r][c])*Real(Ms[b][r][c]) : r, c in [1..2] ])/2;
            G[a][b] := v; G[b][a] := v;
        end for; end for;
        GI := Matrix(Integers(), 4, 4, [ Round(SCL*G[a][b]) : a, b in [1..4] ]);
        for a in [1..4] do for b in [a+1..4] do GI[b][a] := GI[a][b]; end for; end for;
        for a in [1..4] do GI[a][a] +:= 8; end for;
        L := LatticeWithGram(GI);
        cfs := [ (CC ! MoebiusMu(GCD(d, dat`N))) / (d * afs[i][d]) : i in [1..nf] ];
        sub := [ CC | 0 : i in [1..nf] ];
        dsub := [ CC | 0 : i in [1..nf] ];
        P := ShortVectorsProcess(L, Ceiling(SCL*Pb*1.000001));
        while not IsEmpty(P) do
            vv := NextVector(P);
            v := Eltseq(vv);
            if v eq [0,0,0,0] then continue; end if;
            t := &+[ v[a]*v[b]*Qn[a][b] : a, b in [1..4] ];
            if t le 0 then continue; end if;
            n := Integers() ! (d*t);
            if n gt nmax then continue; end if;
            X := v[1]*Xs[1] + v[2]*Xs[2] + v[3]*Xs[3] + v[4]*Xs[4];
            aX := Abs(X); x := 4*pi*aX;
            if x gt cut then continue; end if;
            if Deriv then
                Wv, Wd := WkFastBoth(dat`Wtab, x);
                Xz  := v[1]*Xzs[1]  + v[2]*Xzs[2]  + v[3]*Xzs[3]  + v[4]*Xzs[4];
                Xbz := v[1]*Xbzs[1] + v[2]*Xbzs[2] + v[3]*Xbzs[3] + v[4]*Xbzs[4];
                // d/dz |X| = (conj(X) dX/dz + X d(conj X)/dz) / 2|X|
                du := (ComplexConjugate(X)*Xz + X*Xbz)/(2*aX);
                q := X/aX;
                dq := (Xz - q*du)/aX;
                wv := q^k * Wv;
                dwv := k*q^(k-1)*dq*Wv + q^k*Wd*4*pi*du;
                for i in [1..nf] do
                    sub[i] +:= afs[i][n]*wv;
                    dsub[i] +:= afs[i][n]*dwv;
                end for;
            else
                wv := (X/aX)^k * WkFast(dat`Wtab, x);
                for i in [1..nf] do sub[i] +:= afs[i][n]*wv; end for;
            end if;
            nterms +:= 1;
        end while;
        for i in [1..nf] do tot[i] +:= 2*cfs[i]*sub[i]; end for;
        if Deriv then
            for i in [1..nf] do dtot[i] +:= 2*cfs[i]*dsub[i]; end for;
        end if;
    end for;
    // the derivative comes LAST, so that the existing two-value call sites are unaffected
    return tot, nterms, dtot;
end function;

/////////////////////////////////////////////////////////////////////////////////////////////////
// The W-fixed basis of a quotient

NelBasis := recformat< D, N, W, g, k, prec, cut, cutmax, CC, RR, A, ON, O1, mrep, obas, Wtab,
                       nmax, levels, dat, formsAt, colsAt, wp, labels, halls, alsig, rejsig,
                       newcount, oldcount, volume, z1 >;

// predicted eigenvalue of w_m on the column (level M, coefficients cs, primes T, character chi)
function PredictedEigenvalue(D, M, cs, T, chi, m)
    e := Universe(cs) ! 1;
    for p in PrimeDivisors(m) do
        if D mod p eq 0 then e *:= cs[p];
        elif M mod p eq 0 then e *:= -cs[p];
        else e *:= chi[Index(T, p)]; end if;
    end for;
    return e;
end function;

intrinsic NelsonBasis(D::RngIntElt, N::RngIntElt, W::SetEnum :
                      Precision := 30, Cutoff := 40, Oldforms := true, OrderBasis := [],
                      CarryRejected := false, CoefficientBound := 0, AllowCoarsePrecision := false) -> Rec
{The data needed to evaluate a basis of the W-fixed weight-2 forms of X_0^D(N)/W by Theorem 3.1
 (Nelson, Math. Comp. 84 (2015)) at the given working precision, for cutoffs up to Cutoff: the
 Eichler order (rebuilt from OrderBasis when given, see NelsonEichlerOrderBasis), the real
 splitting, the lattices R^(d) at every level M | N that contributes, the Fourier coefficients of
 the newforms of level DM fixed by W (for M < N through the characters chi of the oldform
 extension; Oldforms := false keeps M = N only), and the table of W_2.  The number of columns must
 equal the genus, which is what checks the Atkin-Lehner eigenvalue convention (a_p at p | D, -a_p
 at p | M, chi(p) at p | N/M); the intrinsic fails otherwise.  With CarryRejected the characters
 that are NOT fixed by W are kept as extra columns, so that NelsonCheckWInvariance can measure
 their predicted eigenvalues (which are not 1) as well.  A request whose evaluations would sum
 more than about a million terms at fewer than 50 digits is refused, since the rounding error of
 the sum then overtakes its truncation error and a residual stalls; AllowCoarsePrecision runs it
 anyway.  Numerical throughout.}
    require IsSquarefree(D*N) and GCD(D, N) eq 1 : "D*N squarefree and coprime is assumed";
    require 1 in W : "W must contain 1";
    k := 2;
    // One evaluation sums about 10 (cut/4pi)^4 sigma_2(DN)/DN terms.  Past about a million of
    // them the rounding error of the sum overtakes its truncation error, and a residual then
    // stops falling with the cutoff while looking like a genuine small singular value, which is
    // the trap this refuses; AllowCoarsePrecision runs it anyway.
    terms := 10 * (Cutoff/(4*Pi(RealField(20))))^4 * (&+[ dv^2 : dv in Divisors(D*N) ])/(D*N);
    require AllowCoarsePrecision or terms lt 10^6 or Precision ge 50 :
        Sprintf("one evaluation sums about %o terms at cutoff %o, which needs 50 digits or more; "
                * "pass Precision := 50, or AllowCoarsePrecision to run at %o anyway",
                Round(terms), Cutoff, Precision);
    PREC := Precision; cutv := Cutoff;
    CC := ComplexField(PREC); RR := RealField(PREC); Ic := CC.1;
    g := GenusShimuraCurveQuotient(D, N, W);
    DN := D*N;
    // Fourier coefficients to the largest d nrd(alpha) the cutoff can reach, plus a margin;
    // CoefficientBound overrides the bound (for measuring that nothing is dropped)
    nmax := Max(60, Ceiling(DN*(cutv/(4*Pi(RR)))^2) + 20);
    if CoefficientBound gt 0 then nmax := CoefficientBound; end if;
    A := QuaternionAlgebra(D);
    // the order is rebuilt inside THIS algebra: a second QuaternionAlgebra(D) is a different
    // structure with no coercion into this one
    if #OrderBasis eq 0 then
        ON := QuaternionOrder(A, N);
    else
        ON := QuaternionOrder([ A ! v : v in OrderBasis ]);
        require forall{ v : v in OrderBasis | IsCoercible(ON, A ! v) } : "OrderBasis does not lie in the order it generates";
    end if;
    require Discriminant(ON) eq DN : "the order does not have discriminant D*N";
    obas := [ Eltseq(A ! x) : x in ZBasis(ON) ];
    O1 := MaximalOrder(ON);
    mrep := FuchsianMatrixRepresentation(A : Precision := PREC);

    // the forms, by level: one entry per newform kept, <M, cs, T, chis, rej>
    Sel := [* *]; newcnt := 0; oldcnt := 0;
    for M in (Oldforms select Divisors(N) else [N]) do
        T := PrimeDivisors(N div M);
        NFs := Newforms(CuspForms(Gamma0(D*M), k));
        for orb in NFs do
            KF := BaseRing(Parent(orb[1]));
            if Level(Parent(orb[1])) ne D*M then continue; end if;
            cj := (Degree(KF) gt 1) select Conjugates(KF.1 : Precision := PREC+10) else [];
            css := [];
            for idx in [1..#orb] do
                f := orb[idx];
                qq := qExpansion(f, nmax+1);
                if Degree(KF) eq 1 then
                    cs := [ CC ! (Rationals() ! Coefficient(qq, n)) : n in [1..nmax] ];
                else
                    cs := [ CC | Evaluate(Polynomial(Rationals(), Eltseq(Coefficient(qq, n))), cj[idx])
                            : n in [1..nmax] ];
                end if;
                Append(~css, cs);
            end for;
            // the Galois conjugates of an orbit are ordered by their coefficients, not by the
            // embedding Magma lists first: the generator of the coefficient field, and with it the
            // order of its conjugates, is not the same from one process to the next, and columns
            // computed in separate processes must agree on which conjugate is which
            Sort(~css, func< a, b | Sign(&+[ Sign(Round(10^12*(Real(a[n]) - Real(b[n])))) * 2^(-n) : n in [1..Min(60, #a)] ]) >);
            for cs in css do
                chis := []; rej := [];
                for bits in CartesianPower([1, -1], #T) do
                    chi := [ bits[i] : i in [1..#T] ];
                    ok := true;
                    for m in W do
                        if m eq 1 then continue; end if;
                        if Abs(PredictedEigenvalue(D, M, cs, T, chi, m) - 1) gt 0.01 then ok := false; break; end if;
                    end for;
                    if ok then Append(~chis, chi); else Append(~rej, chi); end if;
                end for;
                if not CarryRejected then rej := []; end if;
                if #chis eq 0 and #rej eq 0 then continue; end if;
                Append(~Sel, <M, cs, T, chis, rej>);
                if M eq N then newcnt +:= #chis; else oldcnt +:= #chis; end if;
            end for;
        end for;
    end for;
    vprintf NelsonShimizu, 1: "X_0^%o(%o)/%o: genus %o; W-fixed differentials %o from newforms of level %o, %o from oldclasses\n",
        D, N, Sort(SetToSequence(W)), g, newcnt, DN, oldcnt;
    require newcnt + oldcnt eq g :
        Sprintf("the W-fixed differentials found (%o new, %o old) do not number the genus %o%o",
                newcnt, oldcnt, g, Oldforms select "" else "; the oldform part may be needed");

    Wtab := WkTable(k, cutv + 2, RR);
    levels := Sort(SetToSequence({ Integers() | s[1] : s in Sel }));
    dat := AssociativeArray(); formsAt := AssociativeArray(); colsAt := AssociativeArray();
    for M in levels do
        OM := (M eq N) select ON else SuperOrder(ON, O1, M);
        require Discriminant(OM) eq D*M and forall{ x : x in Basis(ON) | A ! x in OM } :
            "the level-M order is not the Eichler order between O_N and the maximal order";
        dat[M] := LevelSetup(D, M, k, PREC, cutv + 2, Wtab, OM, mrep);
        formsAt[M] := [ s[2] : s in Sel | s[1] eq M ];
        cols := []; i := 0;
        for s in Sel do
            if s[1] ne M then continue; end if;
            i +:= 1;
            for chi in s[4] do Append(~cols, <i, chi, true>); end for;
        end for;
        i := 0;
        for s in Sel do
            if s[1] ne M then continue; end if;
            i +:= 1;
            for chi in s[5] do Append(~cols, <i, chi, false>); end for;
        end for;
        colsAt[M] := cols;
    end for;
    wp := AssociativeArray();
    for p in PrimeDivisors(N) do
        if N eq Max(levels) and #levels eq 1 then break; end if;   // no oldclass needs a w_p
        w := ALElement(ON, p, mrep);
        require Norm(w) eq p and (A ! (w*w))/p in ON : "w_p^2 is not p times an element of the order";
        wp[p] := w;
    end for;
    // the Atkin-Lehner signature of every W-fixed column at every Hall divisor, so that a later
    // stage can tell which column belongs to which quotient without re-deriving the convention
    halls := [ m : m in Divisors(DN) | GCD(m, DN div m) eq 1 ];
    labels := []; alsig := []; rejsig := [];
    for M in levels do
        fs := [ s : s in Sel | s[1] eq M ];
        for c in colsAt[M] do
            s := fs[c[1]];
            row := [ Round(Real(PredictedEigenvalue(D, M, s[2], s[3], c[2], m))) : m in halls ];
            if c[3] then
                Append(~labels, Sprintf("level %o form %o chi %o", D*M, c[1], c[2]));
                Append(~alsig, row);
            else
                Append(~rejsig, row);
            end if;
        end for;
    end for;
    vol := RR ! ArithmeticVolume(FuchsianGroup(ON));
    return rec< NelBasis | D := D, N := N, W := W, g := g, k := k, prec := PREC, cut := cutv,
                           cutmax := cutv + 2, CC := CC, RR := RR, A := A, ON := ON, O1 := O1,
                           mrep := mrep, obas := obas, Wtab := Wtab, nmax := nmax, levels := levels,
                           dat := dat, formsAt := formsAt, colsAt := colsAt, wp := wp,
                           labels := labels, halls := halls, alsig := alsig, rejsig := rejsig,
                           newcount := newcnt, oldcount := oldcnt, volume := vol,
                           z1 := CC ! (0.1103 + 1.0411*Ic) >;
end intrinsic;

intrinsic NelsonValues(B::Rec, z1::FldComElt, z2::FldComElt : Cutoff := 0) -> SeqEnum, SeqEnum, RngIntElt
{The values y1 y2 conj(f_i(z1)) f_i(z2) of the W-fixed basis of B (one entry per column, in the
 order of B`labels) by the truncated sum of Theorem 3.1 with the given cutoff (default B`cut);
 then the same for the carried rejected columns (empty unless B was built with CarryRejected),
 and the number of lattice terms summed.  Each column is normalised by its own conj(f_i(z1)), so
 ratios of columns are ratios of forms only across one fixed z1.  Numerical.}
    CC := B`CC; RR := B`RR;
    cut := Cutoff eq 0 select B`cut else Cutoff;
    require cut le B`cutmax : "cutoff above the one the basis was built for";
    A := B`A; mrep := B`mrep; N := B`N;
    z1 := CC ! z1; z := CC ! z2;
    fixed := [ CC | ]; rejv := [ CC | ]; nt := 0;
    for M in B`levels do
        T := PrimeDivisors(N div M);
        subs := Subsets({ p : p in T });
        vS := AssociativeArray();
        for S in subs do
            w := &*([A ! 1] cat [ A ! B`wp[p] : p in S ]);
            Mw := ChangeRing(mrep(w), CC);
            jj := Mw[2][1]*z + Mw[2][2];
            zz := (Mw[1][1]*z + Mw[1][2])/jj;
            vals, n1 := EvalMulti(B`dat[M], z1, zz, B`formsAt[M], RR ! cut, CC);
            nt +:= n1;
            ph := ComplexConjugate(jj)/jj;
            vS[S] := [ ph * v : v in vals ];
        end for;
        for c in B`colsAt[M] do
            tot := CC ! 0;
            for S in subs do
                sg := &*([1] cat [ c[2][Index(T, p)] : p in S ]);
                tot +:= sg * vS[S][c[1]];
            end for;
            if c[3] then Append(~fixed, tot); else Append(~rejv, tot); end if;
        end for;
    end for;
    return fixed, rejv, nt;
end intrinsic;

intrinsic NelsonKernelDerivativeCheck(k::RngIntElt, prec::RngIntElt, cutmax::RngIntElt,
                                      xs::SeqEnum) -> SeqEnum
{For each x in xs, the kernel W_k(x) read off the Chebyshev table, its derivative from the
 differentiated series, and the derivative from a central difference of WkDirect, as a triple.
 This isolates the Chebyshev differentiation from the lattice sum: a convention error in it shows
 up here and nowhere else.}
    RR := RealField(prec);
    Wtab := WkTable(k, cutmax, RR);
    out := [];
    h := RR ! 10^(-(prec div 3));
    for x0 in xs do
        x := RR ! x0;
        W, Wd := WkFastBoth(Wtab, x);
        fd := (WkDirect(k, x+h, RR) - WkDirect(k, x-h, RR))/(2*h);
        Append(~out, <W, WkDirect(k, x, RR), Wd, fd>);
    end for;
    return out;
end intrinsic;

intrinsic NelsonValuesD(B::Rec, z1::FldComElt, z2::FldComElt : Cutoff := 0)
    -> SeqEnum, SeqEnum, RngIntElt
{The values of the W-fixed basis of B at z2, as NelsonValues returns them, together with the
 derivative of each with respect to z2, and the number of lattice terms summed.  The derivative is
 taken inside the lattice sum, so it shares the enumeration of short vectors and costs a fraction
 of a second evaluation, where a difference quotient or a Cauchy circle would cost several.

 WHAT THE SECOND LIST IS.  NelsonValues returns g_i(z) = y1 y2 conj(f_i(z1)) f_i(z2), which is not
 holomorphic because of the factor y2 = Im(z2).  This returns g_i and d g_i/d z2.  The holomorphic
 object is h_i = g_i / Im(z2), and its derivative is recovered as

     h_i'  =  ( d g_i/d z2  +  (i/2) h_i ) / Im(z2),

 since d(Im z)/dz = -i/2.  NelsonDerivedRatio does this and forms the ratio a hyperelliptic
 coordinate needs.  Numerical.}
    CC := B`CC; RR := B`RR;
    cut := Cutoff eq 0 select B`cut else Cutoff;
    require cut le B`cutmax : "cutoff above the one the basis was built for";
    A := B`A; mrep := B`mrep; N := B`N;
    z1 := CC ! z1; z := CC ! z2;
    fixed := [ CC | ]; dfixed := [ CC | ]; nt := 0;
    for M in B`levels do
        T := PrimeDivisors(N div M);
        subs := Subsets({ p : p in T });
        vS := AssociativeArray(); dvS := AssociativeArray();
        for S in subs do
            w := &*([A ! 1] cat [ A ! B`wp[p] : p in S ]);
            Mw := ChangeRing(mrep(w), CC);
            jj := Mw[2][1]*z + Mw[2][2];
            zz := (Mw[1][1]*z + Mw[1][2])/jj;
            vals, n1, dvals := EvalMulti(B`dat[M], z1, zz, B`formsAt[M], RR ! cut, CC : Deriv);
            nt +:= n1;
            ph := ComplexConjugate(jj)/jj;
            // Mw is a real matrix, so conj(jj) is a function of conj(z) alone and drops out of
            // d/dz; the two surviving pieces are d/dz of 1/jj and the inner substitution
            dph := -Mw[2][1]*ComplexConjugate(jj)/jj^2;
            dzz := (Mw[1][1]*Mw[2][2] - Mw[1][2]*Mw[2][1])/jj^2;
            vS[S] := [ ph * v : v in vals ];
            dvS[S] := [ dph*vals[i] + ph*dvals[i]*dzz : i in [1..#vals] ];
        end for;
        for c in B`colsAt[M] do
            tot := CC ! 0; dtot := CC ! 0;
            for S in subs do
                sg := &*([1] cat [ c[2][Index(T, p)] : p in S ]);
                tot +:= sg * vS[S][c[1]];
                dtot +:= sg * dvS[S][c[1]];
            end for;
            if c[3] then Append(~fixed, tot); Append(~dfixed, dtot); end if;
        end for;
    end for;
    return fixed, dfixed, nt;
end intrinsic;

intrinsic NelsonDerivedRatio(B::Rec, z1::FldComElt, z2::FldComElt : Cutoff := 0)
    -> SeqEnum, SeqEnum, FldComElt, FldComElt
{The forms h_i = g_i / Im(z2) of the W-fixed basis at z2, their derivatives h_i', and for a
 genus-two quotient the hyperelliptic coordinate x = h_2/h_1 together with dx/dz2.  This is the
 pair a hyperelliptic equation needs: with omega_i = h_i dz, the curve y^2 = F(x) has
 x = h_2/h_1 and y = (dx/dz)/h_1, so F is then a linear system in its seven coefficients.}
    CC := B`CC;
    g, dg, _ := NelsonValuesD(B, z1, z2 : Cutoff := Cutoff);
    y := Imaginary(CC ! z2);
    h := [ CC | t/y : t in g ];
    hd := [ CC | (dg[i] + (CC.1/2)*h[i])/y : i in [1..#g] ];
    require #h ge 2 : "a hyperelliptic coordinate needs at least two forms";
    x := h[2]/h[1];
    dx := (hd[2]*h[1] - h[2]*hd[1])/h[1]^2;
    return h, hd, x, dx;
end intrinsic;

intrinsic NelsonModelFromDifferentials(D::RngIntElt, N::RngIntElt, W::SetEnum :
            Precision := 60, Cutoff := 46, Samples := 0, Discriminants := [],
            DiscriminantBound := 400, PrimeCap := 300, OrderBasis := [],
            AllowCoarsePrecision := false) -> CrvHyp, Rec
{A hyperelliptic model over Q of a GENUS-TWO Atkin-Lehner quotient, from its differentials alone,
 with the diagnostics record second.  Unlike NelsonHyperellipticModel this needs no marked
 Atkin-Lehner involution, so it reaches the quotients where W is already the full star.

 In genus two the canonical map IS the hyperelliptic map: with omega_i = h_i dz the coordinate
 x = h_2/h_1 and y = (dx/dz)/h_1 satisfy y^2 = F(x), which is LINEAR in the seven coefficients of
 F, so seven points determine it and the others test it.  The derivative comes from inside the
 lattice sum (NelsonDerivedRatio), one evaluation per point.

 The coefficients are complex, the basis being a complex one, so the coordinate is then fixed by
 sending three CM points defined over Q, in three DISTINCT fibres of the hyperelliptic map, to 0,
 1 and infinity (rational CM points pair up under the hyperelliptic involution, and the pairing is
 returned, being the degree-two divisor an exact construction would prescribe); in that coordinate every rational
 point has rational coordinate, so the sextic has rational coefficient ratios, of small height.
 Recognising the absolute invariants instead does not work: their heights are far beyond what the
 values carry.  The discriminants are found here rather than supplied, since normalising on a point
 that is NOT rational yields a plausible wrong curve that nothing downstream detects; pass
 Discriminants to override, and DiscriminantBound to widen the search.

 What is left over is one scalar: the values carry an unknown constant of the base point,
 h_i = K_i f_i, which rescales x harmlessly but multiplies y by K_2/K_1^2, and normalising by a
 coefficient discards it.  The quadratic twist is therefore read off the point counts against the
 newform orbit, matched on |Tr a_p| so that no Atkin-Lehner sign convention enters.  Those counts
 fix the twist and confirm the isogeny class; they could NOT select the curve, since every
 principal polarization on the surface shares an L-function.

 Numerical.  Fails with a message at the first step that does not meet its own test: the points
 held back from the fit, the collinearity of the rational CM points, the rationality of the
 coefficients, or the counts.}
    require IsSquarefree(D*N) and GCD(D, N) eq 1 : "D*N squarefree and coprime is assumed";
    nc := (Samples eq 0) select 10 else Samples;
    require nc ge 9 : "at least nine points, so that two are held back";

    // which discriminants give a CM point defined over Q, and which of those sit at an
    // Atkin-Lehner fixed point, where x = h_2/h_1 is 0/0
    rds := Discriminants;
    X := CreateShimuraQuot(D, N, W);
    X`g := GenusShimuraCurveQuotient(D, N, W);
    require X`g eq 2 : Sprintf("this is for genus two; the quotient has genus %o", X`g);
    if IsEmpty(rds) then
        for d in [-3 .. -DiscriminantBound by -1] do
            if d mod 4 notin {0, 1} then continue; end if;
            dg := 0;
            try dg := DegreeOfFieldOfDefinitionOfCMPoint(X, d); catch e continue; end try;
            if dg eq 1 then Append(~rds, d); end if;
        end for;
    end if;
    alfix := func< m | m eq 1 select {} else
                       (m eq 2 select {-4, -8} else
                       ((m mod 4) eq 3 select {-m, -4*m} else {-4*m})) >;
    isbr := [ exists{ m : m in W | d in alfix(m) } : d in rds ];
    // the three that fix the coordinate must be ordinary and have distinct coordinates, so the
    // ordinary ones are tried first
    rds := [ rds[i] : i in [1..#rds] | not isbr[i] ] cat [ rds[i] : i in [1..#rds] | isbr[i] ];
    isbr := [ exists{ m : m in W | d in alfix(m) } : d in rds ];
    require #[ i : i in [1..#rds] | not isbr[i] ] ge 3 :
        Sprintf("only %o CM points of X are both defined over Q and away from an Atkin-Lehner "
                cat "fixed point; three are needed to fix the coordinate", 
                #[ i : i in [1..#rds] | not isbr[i] ]);

    B := NelsonBasis(D, N, W : Precision := Precision, Cutoff := Cutoff, OrderBasis := OrderBasis,
                               AllowCoarsePrecision := AllowCoarsePrecision);
    require B`g eq 2 : "the basis does not have genus two";
    CC := B`CC; RR := RealField(Precision); A := B`A; O := B`ON;
    PQ := PolynomialRing(Rationals()); T := PQ.1;
    diag := rec< recformat< discs, isbranch, cmx, fibres, normalisers, spread, worst, tol, bound,
                            margins, sextic, invariants, twist, counts, orbit > | >;

    bas := [ A ! t : t in Basis(O) ];
    Lsp := Kernel(Transpose(Matrix(Integers(), [[ Integers() ! Trace(t) : t in bas ]])));
    bl := [ &+[ b[i]*bas[i] : i in [1..4] ] : b in Basis(Lsp) ];
    Qf := Matrix([[ Norm(t+u) - Norm(t) - Norm(u) : u in bl ] : t in bl ]);
    function CMPoint(d)
        found, v := FindLambda(Qf, -d, O, bl);
        if not found then return false, CC ! 0; end if;
        elt := &+[ v[i]/2*bl[i] : i in [1..#bl] ] + (((-d) mod 4 eq 3) select A ! (1/2) else A ! 0);
        th := A ! (2*elt - Trace(elt));
        if th^2 ne d then return false, CC ! 0; end if;
        Mm := ChangeRing(B`mrep(th), CC);
        a := Mm[1][1]; bb := Mm[1][2]; c := Mm[2][1];
        rt := Sqrt(CC ! (a^2 + bb*c));
        z := (a + rt)/c;
        if Imaginary(z) lt 0 then z := (a - rt)/c; end if;
        return true, z;
    end function;

    cmx := [ CC | ]; cmy := [ CC | ]; keep := []; kbr := [];
    for i in [1..#rds] do
        ok, z := CMPoint(rds[i]);
        if not ok then continue; end if;
        h, hd, xv, dxv := NelsonDerivedRatio(B, B`z1, z : Cutoff := Cutoff);
        if isbr[i] then
            // both forms vanish; an Atkin-Lehner involution has order two, so each vanishes to
            // order exactly one and the coordinate is the limit h_2'/h_1'
            Append(~cmx, hd[2]/hd[1]); Append(~cmy, CC ! 0);
        else
            Append(~cmx, xv); Append(~cmy, dxv/h[1]);
        end if;
        Append(~keep, rds[i]); Append(~kbr, isbr[i]);
    end for;
    require #[ i : i in [1..#keep] | not kbr[i] ] ge 3 :
        "fewer than three ordinary rational CM points embed in this order";
    diag`discs := keep; diag`isbranch := kbr; diag`cmx := cmx;
    // a rational point must lie on the image of the real locus, which is a line through the origin
    // in this normalisation; this is NECESSARY, not sufficient
    rsp := [ Imaginary(cmx[i])/Real(cmx[i]) : i in [1..#cmx] ];
    diag`spread := Maximum(rsp) - Minimum(rsp);

    // the sextic.  Nodes are chosen FARTHEST APART and the ordinary CM points join the sample
    // points: the sample coordinates cluster, and a degree-six Vandermonde on close nodes discards
    // most of the digits the evaluation delivered.
    pts := NelsonSamplePoints(B, nc);
    alx := [ CC | ]; aly := [ CC | ];
    for j in [1..nc] do
        h, hd, xv, dxv := NelsonDerivedRatio(B, B`z1, pts[j] : Cutoff := Cutoff);
        Append(~alx, xv); Append(~aly, dxv/h[1]);
    end for;
    for i in [1..#keep] do
        if not kbr[i] then Append(~alx, cmx[i]); Append(~aly, cmy[i]); end if;
    end for;
    nu := 7;
    sel := [ 1 ];
    for t in [2..nu] do
        bd := RR ! -1; bi := 0;
        for i in [1..#alx] do
            if i in sel then continue; end if;
            dd := Minimum([ Abs(alx[i] - alx[j]) : j in sel ]);
            if dd gt bd then bd := dd; bi := i; end if;
        end for;
        Append(~sel, bi);
    end for;
    held := [ i : i in [1..#alx] | i notin sel ];
    Mv := Matrix(CC, nu, nu, [ alx[sel[i]]^(k-1) : k in [1..nu], i in [1..nu] ]);
    require Abs(Determinant(Mv)) ge 10^(-Precision div 2) : "the chosen nodes do not determine F";
    cs := Eltseq(Mv^-1 * Matrix(CC, nu, 1, [ aly[sel[i]]^2 : i in [1..nu] ]));
    PCx := PolynomialRing(CC); xx := PCx.1;
    F := &+[ cs[k]*xx^(k-1) : k in [1..nu] ];
    // the residual is scaled by the size of the TERMS, not by |y|^2: a held-back point near a
    // Weierstrass point has y near zero and would report a huge residual on a perfect fit
    worst := RR ! 0;
    for i in held do
        sc := &+[ Abs(cs[k])*Abs(alx[i])^(k-1) : k in [1..nu] ] + Abs(aly[i])^2;
        worst := Max(worst, Abs(Evaluate(F, alx[i]) - aly[i]^2)/sc);
    end for;
    diag`worst := worst;

    // The three points that fix the coordinate must lie in DISTINCT fibres of the hyperelliptic
    // map.  Rational CM points pair up under the hyperelliptic involution -- at 2002 the points of
    // discriminant -148 and -232 share an x, and so do -163 and -312 -- and two of the three
    // sharing an x gives a Mobius map of enormous scale, after which every coefficient but one
    // recognises as zero.  So the three are chosen farthest apart, and the pairing is recorded,
    // since it is exactly the degree-two divisor the exact (Borcherds) route needs to prescribe.
    ord := [ i : i in [1..#keep] | not kbr[i] ];
    sep := RR ! 10^(-8);
    fibres := [];
    for i in ord do
        placed := false;
        for f in [1..#fibres] do
            if Abs(cmx[i] - cmx[fibres[f][1]]) le sep*(1 + Abs(cmx[i])) then
                Append(~fibres[f], i); placed := true; break;
            end if;
        end for;
        if not placed then Append(~fibres, [i]); end if;
    end for;
    diag`fibres := [ [ keep[i] : i in f ] : f in fibres ];
    require #fibres ge 3 :
        Sprintf("the rational CM points away from Atkin-Lehner fixed points lie in only %o "
                cat "fibre(s) of the hyperelliptic map (%o); three distinct fibres are needed to "
                cat "fix the coordinate", #fibres, diag`fibres);
    reps := [ f[1] : f in fibres ];
    pick := [ reps[1] ];
    for t in [2..3] do
        bd := RR ! -1; bi := 0;
        for i in reps do
            if i in pick then continue; end if;
            dd := Minimum([ Abs(cmx[i] - cmx[j]) : j in pick ]);
            if dd gt bd then bd := dd; bi := i; end if;
        end for;
        Append(~pick, bi);
    end for;
    x1 := cmx[pick[1]]; x2 := cmx[pick[2]]; x3 := cmx[pick[3]];
    diag`normalisers := [ keep[i] : i in pick ];
    Ac := x2 - x3; Bc := -x1*(x2 - x3); Cc := x2 - x1; Dc := -x3*(x2 - x1);
    PCu := PolynomialRing(CC); uu := PCu.1;
    G := &+[ cs[k+1] * (-Dc*uu + Bc)^k * (Cc*uu - Ac)^(6-k) : k in [0..6] ];
    gc := [ Coefficient(G, k) : k in [0..6] ];
    _, imax := Maximum([ Abs(t) : t in gc ]);
    gn := [ t/gc[imax] : t in gc ];
    // BOTH the tolerance and the denominator bound come from the accuracy the fit MEASURED at the
    // points held back, never from the working precision: BestApproximation(x, Q) returns a
    // fraction within 1/Q^2 of any real whatever, so a bound taken from the precision guarantees a
    // spurious fit.  Q is chosen so that a chance fit would be four orders of magnitude worse.
    tol := Maximum(100*worst, RR ! 10^(-Precision + 5));
    Qb := Maximum(10, Floor(Sqrt(1/tol)/100));
    diag`tol := tol; diag`bound := Qb;
    qs := []; mg := [];
    for k in [1..nu] do
        z := gn[k];
        q := BestApproximation(Real(z), Qb);
        err := Max(Abs(Imaginary(z)), Abs(RR ! q - Real(z)))/(1 + Abs(z));
        require err le tol :
            Sprintf("coefficient %o is not rational at the accuracy achieved (%o against a "
                    cat "tolerance of %o): raise the cutoff AND the precision", k-1, err, tol);
        Append(~qs, q); Append(~mg, (err gt 0) select (RR ! 1/Qb^2)/err else RR ! 0);
    end for;
    diag`margins := mg;
    den := LCM([ Denominator(q) : q in qs ]);
    GQ := &+[ (Integers() ! (den*qs[k+1])) * T^k : k in [0..6] ];
    diag`sextic := GQ;
    require Discriminant(GQ) ne 0 : "the recognised sextic is singular";
    diag`invariants := G2Invariants(HyperellipticCurve(GQ));

    // the twist, from the counts against the newform orbit
    supp := Set(PrimeDivisors(D*N))
            join Set(PrimeDivisors(Integers() ! (den^12 * Discriminant(GQ))));
    for c in Coefficients(GQ) do
        if c ne 0 then supp join:= Set(PrimeDivisors(Numerator(c))); end if;
        supp join:= Set(PrimeDivisors(Denominator(c)));
    end for;
    supp := Sort(SetToSequence(supp));
    function CountAt(poly, p)
        dd := LCM([ Denominator(c) : c in Coefficients(poly) ] cat [1]);
        if dd mod p eq 0 then return -1; end if;
        gg := ChangeRing(PQ ! (dd^2 * poly), GF(p));
        if Degree(gg) lt Degree(PQ ! (dd^2*poly)) - 1 or Degree(gg) lt 3
           or Discriminant(gg) eq 0 then return -1; end if;
        return #Points(HyperellipticCurve(gg));
    end function;
    cnt := AssociativeArray();
    for p in PrimesInInterval(3, PrimeCap) do
        if (D*N) mod p ne 0 then
            c := CountAt(GQ, p);
            if c ge 0 then cnt[p] := c; end if;
        end if;
    end for;
    gp := Sort([ p : p in Keys(cnt) ]);
    hit := 0; traces := AssociativeArray();
    for Fm in NewformDecomposition(NewSubspace(CuspidalSubspace(ModularSymbols(D*N, 2)))) do
        if Dimension(Fm) ne 4 then continue; end if;
        tr := AssociativeArray(); good := true;
        for p in gp do
            // Trace(HeckeOperator) is TWICE Tr(a_p): T_p acts on the four-dimensional homology
            // as a_p on K^2.  The orbit is matched on |Tr a_p|, which is convention-free.
            t := (Integers() ! Trace(HeckeOperator(Fm, p))) div 2;
            tr[p] := t;
            if Abs(p + 1 - cnt[p]) ne Abs(t) then good := false; break; end if;
        end for;
        if good then hit := Fm; traces := tr; break; end if;
    end for;
    require hit cmpne 0 :
        "no dimension-two newform orbit of this level matches |Tr a_p| at every good prime, so the "
        cat "curve is not in the isogeny class of any of them";
    diag`orbit := hit;
    rows := []; rhs := []; used := [];
    for p in gp do
        t := traces[p];
        if t eq 0 then continue; end if;         // carries no information about the character
        b := (p + 1 - cnt[p] eq t) select 0 else 1;
        Append(~rows, [ GF(2) | (LegendreSymbol(q, p) eq 1 select 0 else 1) : q in supp ]
                      cat [ GF(2) | (p mod 4 eq 1 select 0 else 1) ]);
        Append(~rhs, GF(2) ! b); Append(~used, p);
    end for;
    require #used ge #supp + 4 : "too few informative primes for the twist; raise PrimeCap";
    Mt := Matrix(GF(2), #rows, #supp+1, &cat rows);
    ok2, vv := IsConsistent(Transpose(Mt), Vector(GF(2), rhs));
    require ok2 : "the Legendre symbols are inconsistent: no single twist fits the counts";
    require Rank(Mt) eq #supp+1 : "the twist is not determined uniquely; raise PrimeCap";
    sgn := (vv[#supp+1] eq GF(2)!1) select -1 else 1;
    tw := sgn * (&*([Integers() | 1] cat [ supp[i] : i in [1..#supp] | vv[i] eq GF(2)!1 ]));
    diag`twist := tw;
    Cm := ReducedMinimalWeierstrassModel(IntegralModel(HyperellipticCurve(tw*GQ)));
    bad := [ p : p in gp | CountAt(tw*GQ, p) ne p + 1 - traces[p] ];
    diag`counts := <#gp, #bad>;
    require IsEmpty(bad) :
        Sprintf("the twisted model misses the newform's count at %o of %o good primes", #bad, #gp);
    return Cm, diag;
end intrinsic;

intrinsic NelsonChooseBasePoint(~B::Rec : Candidates := 6, Cutoff := 0)
{Replace the base point z1 of B (the point at which every column is normalised through
 conj(f_i(z1))) by the best of Candidates points: the default z1 and sample points, scored by the
 smallest diagonal value y^2 |f_i(z)|^2 over the columns relative to its largest value among
 the candidates.  At a base point where a diagonal value is small -- a point near a zero of one of
 the forms -- every value of that column is a small difference of large terms, so its relative
 truncation error is far larger than elsewhere at the same cutoff; which points these are depends
 on the Eichler order.  Costs Candidates evaluations.  Numerical.}
    cands := [ B`z1 ] cat NelsonSamplePoints(B, Candidates - 1);
    diag := [];
    for z in cands do
        v := NelsonValues(B, z, z : Cutoff := Cutoff);
        Append(~diag, [ Abs(x) : x in v ]);
    end for;
    best := 1; bestv := B`RR ! -1;
    for j in [1..#cands] do
        m := Min([ diag[j][i] / Max([ diag[l][i] : l in [1..#cands] ]) : i in [1..B`g] ]);
        if m gt bestv then bestv := m; best := j; end if;
    end for;
    vprintf NelsonShimizu, 1: "base point %o of %o chosen, smallest relative diagonal %o\n",
        best, #cands, RealField(6) ! bestv;
    B`z1 := cands[best];
end intrinsic;

intrinsic NelsonSamplePoints(B::Rec, n::RngIntElt) -> SeqEnum
{n points of the upper half-plane spread by the square root of the arithmetic volume of the
 group of B.  Clustered points give rows that are nearly proportional, and a monomial matrix of
 effective rank near the number of columns always has a quadric through it, so the spread is what
 makes Noether's count readable.}
    CC := B`CC; RR := B`RR; Ic := CC.1;
    sp := Sqrt(B`volume);
    return [ CC | (sp*0.73*Cos(RR ! (1.1*j^2 + 0.3*j)))
                  + (0.55 + sp*0.42*(1 + Sin(RR ! (0.7*j^2 + 1.9*j))))*Ic : j in [1..n] ];
end intrinsic;

intrinsic NelsonAtkinLehnerElement(B::Rec, m::RngIntElt) -> AlgQuatElt
{An element of the Eichler order of B of reduced norm m normalising it, for a Hall divisor m of DN.}
    require m in B`halls : "m must be a Hall divisor of D*N";
    return ALElement(B`ON, m, B`mrep);
end intrinsic;

intrinsic NelsonCheckWInvariance(B::Rec, pts::SeqEnum : Cutoff := 0) -> FldReElt, SeqEnum
{For each generator m of W and each point z of pts, the ratio of the value of every column at
 w_m z, times conj(j)/j, to its value at z, divided by the predicted eigenvalue (1 for the W-fixed
 columns, the character value for carried rejected ones).  Returns the largest deviation of these
 ratios from 1 and the list of <m, z, ratios>.  A wrong relative scale c between a newform and
 its pullback would show as c and 1/c in place of +1 and -1.  Numerical.}
    CC := B`CC;
    gens := ALGenerators(B`W);
    worst := B`RR ! 0; out := [];
    for m in gens do
        wm := ALElement(B`ON, m, B`mrep);
        Mw := ChangeRing(B`mrep(wm), CC);
        for z0 in pts do
            z := CC ! z0;
            v1, r1 := NelsonValues(B, B`z1, z : Cutoff := Cutoff);
            jj := Mw[2][1]*z + Mw[2][2]; zz := (Mw[1][1]*z + Mw[1][2])/jj;
            v2, r2 := NelsonValues(B, B`z1, zz : Cutoff := Cutoff);
            ph := ComplexConjugate(jj)/jj;
            im := Index(B`halls, m);
            rat := [ v2[c]*ph/v1[c]/B`alsig[c][im] : c in [1..#v1] ]
                   cat [ r2[c]*ph/r1[c]/B`rejsig[c][im] : c in [1..#r1] ];
            for r in rat do worst := Max(worst, Abs(r - 1)); end for;
            Append(~out, <m, z, rat>);
        end for;
    end for;
    return worst, out;
end intrinsic;

/////////////////////////////////////////////////////////////////////////////////////////////////
// Noether's count

intrinsic NoetherQuadricCount(vals::SeqEnum, g::RngIntElt : MinGap := 10^4, ZeroTol := 10^(-4)) -> RngIntElt, SeqEnum, FldReElt
{Given values of a basis of g differentials at several points (vals[j] the g values at the j-th
 point), the number of quadrics through the canonical image, read off the singular values of the
 matrix of quadratic monomials: rows and columns are normalised, the matrix is realified, and the
 count is the number of singular values below ZeroTol, accepted only when that number is one a
 canonical image of genus g can have -- (g-2)(g-3)/2 if the curve is not hyperelliptic,
 (g-1)(g-2)/2 if it is -- and the ratio of the first value above them to the last below them
 reaches MinGap; otherwise the count is -1, inconclusive.  Returns the count, the singular values
 (ascending) and that ratio.  Counting the values below a tolerance rather than looking for the
 largest ratio anywhere is what keeps a cluster of truncation-level values that differ among
 themselves from being split, which at genus 4 would report one quadric for a hyperelliptic curve.
 Reading rule: a small singular value is a small number.  A value that falls as the cutoff rises
 is the truncation of the sum; a value that does NOT fall may be genuine or may be a limit of the
 working precision, and only repeating the count at another precision separates those -- at one
 base of this project the smallest value sat at 1.9e-3 at three cutoffs at 50 digits and at
 2.0e-10 at 80, so a value anywhere in the spectrum can be precision-limited and stable across
 cutoffs.  See paper/nelson-method.tex, "What a numerical kernel dimension proves", and use
 NelsonHyperellipticVerdict's ConfirmPrecision.  Genus 3: one relation means hyperelliptic, none
 means not; genus 4: three against one.  Numerical.}
    npts := #vals;
    CC := Universe(vals[1]); RR := RealField(Precision(CC));
    mons := [ [a,b] : b in [a..g], a in [1..g] ];
    require npts ge #mons : "need at least as many points as quadratic monomials";
    M := Matrix(CC, npts, #mons, &cat[ [ vals[j][m[1]]*vals[j][m[2]] : m in mons ] : j in [1..npts] ]);
    for i in [1..npts] do
        nm := Sqrt(&+[ Abs(M[i][c])^2 : c in [1..#mons] ]);
        for c in [1..#mons] do M[i][c] := M[i][c]/nm; end for;
    end for;
    for c in [1..#mons] do
        nm := Sqrt(&+[ Abs(M[i][c])^2 : i in [1..npts] ]);
        for i in [1..npts] do M[i][c] := M[i][c]/nm; end for;
    end for;
    MR := Matrix(RR, 2*npts, 2*#mons, [ RR ! 0 : i in [1..4*npts*#mons] ]);
    for i in [1..npts] do for j in [1..#mons] do
        a := Real(M[i][j]); b := Imaginary(M[i][j]);
        MR[i][j] := a; MR[i][#mons+j] := -b; MR[npts+i][j] := b; MR[npts+i][#mons+j] := a;
    end for; end for;
    Sig := SingularValueDecomposition(MR);
    r := Min(Nrows(MR), Ncols(MR));
    sv := Reverse([ Abs(Sig[i][i]) : i in [1..r] ]);
    sv := [ sv[i] : i in [1..#sv by 2] ];      // realification doubles every singular value
    // The count is the number of singular values below ZeroTol, and it is accepted only if that
    // number is one a canonical image of genus g can have -- (g-2)(g-3)/2 when the curve is not
    // hyperelliptic, (g-1)(g-2)/2 when it is -- and if the gap above them reaches MinGap.
    // Reading the count off the largest ratio anywhere instead can land inside the cluster of
    // truncation-level values when they differ among themselves, which at genus 4 reports one
    // quadric for a hyperelliptic curve.
    nonhyp := (g-2)*(g-3) div 2; hyp := (g-1)*(g-2) div 2;
    n := #[ i : i in [1..#sv] | sv[i] lt ZeroTol ];
    if n eq 0 then return (nonhyp eq 0 select 0 else -1), sv, sv[1]/ZeroTol; end if;
    if n ge #sv then return -1, sv, RR ! 0; end if;
    // exact input (a curve given by equations rather than by evaluated values) makes the
    // relations vanish exactly, and an exact zero is the cleanest gap there is
    gap := sv[n] eq 0 select RR ! 10^Precision(RR) else sv[n+1]/sv[n];
    if not (n in { nonhyp, hyp }) or gap lt MinGap then return -1, sv, gap; end if;
    return n, sv, gap;
end intrinsic;

intrinsic NelsonHyperellipticVerdict(D::RngIntElt, N::RngIntElt, W::SetEnum :
                                     Precision := 30, Cutoff := 40, Points := 18, Oldforms := true,
                                     CheckW := 0, OrderBasis := [], ChooseBasePoint := true,
                                     AllowCoarsePrecision := false, ConfirmPrecision := 0,
                                     AllowUnconfirmed := false) -> MonStgElt, Rec
{Noether's count of canonical quadrics on X_0^D(N)/W, genus at least 3, with the differentials
 evaluated by Theorem 3.1 at Points sample points: "HYPERELLIPTIC", "NOT hyperelliptic" or
 "INCONCLUSIVE", and a record with the singular values, the gap, the number of relations, the
 lattice-term count, and with CheckW := n the largest deviation of the W-invariance of every
 column measured at n points.  The base point is chosen by NelsonChooseBasePoint unless
 ChooseBasePoint is false.
 ConfirmPrecision := m repeats the count at m digits, at the same points and the same base point,
 and reports whether the smallest singular value and the count agree; when they disagree the
 verdict returned is the one at the higher precision, which should itself be confirmed by a
 further raise.  A NOT hyperelliptic verdict is REFUSED without that confirmation (or without
 AllowUnconfirmed), because it is the reading a precision floor can fake: at one base of this
 project the smallest singular value sat at 1.9e-3 at three cutoffs at 50 digits, agreeing to
 eleven digits, and was 2.0e-10 at 80 digits -- not hyperelliptic became hyperelliptic.  Taking
 the count at two cutoffs is NOT sufficient on its own, since two cutoffs can share one precision
 floor.  A verdict is a reading of numerical singular values, not a proof.}
    B := NelsonBasis(D, N, W : Precision := Precision, Cutoff := Cutoff, Oldforms := Oldforms,
                              OrderBasis := OrderBasis, CarryRejected := CheckW gt 0,
                              AllowCoarsePrecision := AllowCoarsePrecision);
    require B`g ge 3 : "Noether's count needs genus at least 3";
    if ChooseBasePoint then NelsonChooseBasePoint(~B); end if;
    pts := [ B`z1 ] cat [ z : z in NelsonSamplePoints(B, Points) | z ne B`z1 ][1..Points - 1];
    vals := []; nt := 0;
    for z in pts do
        v, _, n := NelsonValues(B, B`z1, z);
        Append(~vals, v); nt +:= n;
    end for;
    cnt, sv, gap := NoetherQuadricCount(vals, B`g);
    g := B`g;
    verdict := cnt eq (g-1)*(g-2) div 2 select "HYPERELLIPTIC"
               else (cnt eq (g-2)*(g-3) div 2 select "NOT hyperelliptic" else "INCONCLUSIVE");
    wdev := B`RR ! 0;
    if CheckW gt 0 then
        wdev := NelsonCheckWInvariance(B, pts[1..Min(CheckW, #pts)]);
    end if;
    vprintf NelsonShimizu, 1: "singular values %o\n   relations %o, gap %o, W-check %o: %o\n",
        [ RealField(6) ! x : x in sv ], cnt, RealField(6) ! gap, RealField(6) ! wdev, verdict;
    confirmed := false; svc := []; cntc := 0; verdictc := "";
    if ConfirmPrecision gt 0 then
        B2 := NelsonBasis(D, N, W : Precision := ConfirmPrecision, Cutoff := Cutoff,
                                   Oldforms := Oldforms, OrderBasis := B`obas,
                                   AllowCoarsePrecision := AllowCoarsePrecision);
        B2`z1 := B2`CC ! B`z1;          // the same base point, so the two counts are comparable
        valsc := [];
        for z in pts do
            v := NelsonValues(B2, B2`z1, B2`CC ! z);
            Append(~valsc, v);
        end for;
        cntc, svc, gapc := NoetherQuadricCount(valsc, g);
        verdictc := cntc eq (g-1)*(g-2) div 2 select "HYPERELLIPTIC"
                    else (cntc eq (g-2)*(g-3) div 2 select "NOT hyperelliptic" else "INCONCLUSIVE");
        moved := svc[1] eq 0 or sv[1] eq 0 select true
                 else Maximum(svc[1]/sv[1], sv[1]/svc[1]) gt 10;
        confirmed := (verdictc eq verdict) and not moved;
        vprintf NelsonShimizu, 1: "at %o digits: smallest %o against %o, %o -- %o\n",
            ConfirmPrecision, RealField(6) ! svc[1], RealField(6) ! sv[1], verdictc,
            confirmed select "confirmed" else "NOT confirmed: the lower precision was inadequate";
        if not confirmed then verdict := verdictc; end if;
    end if;
    require verdict ne "NOT hyperelliptic" or ConfirmPrecision gt 0 or AllowUnconfirmed :
        "a NOT hyperelliptic reading is the one a precision floor can fake, and two cutoffs can "
        * "share one floor: pass ConfirmPrecision := (more digits) to repeat the count there, or "
        * "AllowUnconfirmed to take it as it stands";
    R := rec< recformat< singular, gap, relations, terms, wcheck, basis, values, points,
                         confirmed, singularatconfirm, relationsatconfirm, verdictatconfirm > |
              singular := sv, gap := gap, relations := cnt, terms := nt, wcheck := wdev,
              basis := B, values := vals, points := pts, confirmed := confirmed,
              singularatconfirm := svc, relationsatconfirm := cntc, verdictatconfirm := verdictc >;
    return verdict, R;
end intrinsic;

/////////////////////////////////////////////////////////////////////////////////////////////////
// Values on circles, for the model.  The canonical image of a hyperelliptic curve carries no
// moduli, so the values of the differentials alone cannot give the equation; one derivative does,
// and it is taken by Cauchy's formula on a circle rather than by finite differences, which would
// cost a third of the digits.  See paper/nelson-method.tex, "The model".

NelVals := recformat< D, N, W, g, prec, cut, rho, ncirc, npc, halls, alsig, pts, owner, isctr, vals, obas >;

intrinsic NelsonCircleValues(B::Rec, centres::SeqEnum, np::RngIntElt, rho::FldReElt :
                             Cutoff := 0, Range := [0, 0]) -> Rec
{The values of the basis of B at each centre and at np equally spaced points of the circle of
 radius rho around it, as a record that WriteNelsonValues stores and NelsonHyperellipticModel
 consumes.  Range := [a, b] evaluates only the points with index a..b (the point list is centre 1,
 its circle, centre 2, ...) leaving the others zero, so that chunks computed in separate processes
 against the SAME order (B built with OrderBasis) can be combined with NelsonMergeValues.  Numerical.}
    CC := B`CC; Ic := CC.1;
    nc := #centres;
    pts := [ CC | ]; owner := []; isctr := [];
    for c in [1..nc] do
        Append(~pts, CC ! centres[c]); Append(~owner, c); Append(~isctr, 1);
        for t in [0..np-1] do
            Append(~pts, CC ! centres[c] + rho*Exp(2*Pi(CC)*Ic*t/np));
            Append(~owner, c); Append(~isctr, 0);
        end for;
    end for;
    js := Range[1] eq 0 select 1 else Range[1];
    je := Range[2] eq 0 select #pts else Range[2];
    require js ge 1 and je le #pts and js le je : "Range out of bounds";
    vals := [ [ CC | 0 : i in [1..B`g] ] : j in [1..#pts] ];
    t0 := Cputime();
    for j in [js..je] do
        v := NelsonValues(B, B`z1, pts[j] : Cutoff := Cutoff);
        vals[j] := v;
        if (j - js + 1) mod 10 eq 0 then
            vprintf NelsonShimizu, 1: "  %o/%o  %o s\n", j - js + 1, je - js + 1, RealField(6) ! Cputime(t0);
        end if;
    end for;
    return rec< NelVals | D := B`D, N := B`N, W := Sort(SetToSequence(B`W)), g := B`g, prec := B`prec,
                          cut := Cutoff eq 0 select B`cut else Cutoff, rho := B`RR ! rho,
                          ncirc := nc, npc := np, halls := B`halls, alsig := B`alsig,
                          pts := pts, owner := owner, isctr := isctr, vals := vals, obas := B`obas >;
end intrinsic;

intrinsic NelsonMergeValues(recs::SeqEnum) -> Rec
{Combine value records of the same base, order, precision, cutoff and radius: records with the
 same point list (chunks of one circle set) are merged entry by entry, taking the nonzero value;
 records with different centres are concatenated as further circles.}
    require #recs ge 1 : "nothing to merge";
    V := recs[1];
    for i in [2..#recs] do
        U := recs[i];
        require U`D eq V`D and U`N eq V`N and U`W eq V`W and U`prec eq V`prec and U`cut eq V`cut
                and U`npc eq V`npc and U`obas eq V`obas and U`halls eq V`halls and U`alsig eq V`alsig :
            "records differ in base, order, precision, cutoff or circle size";
        if U`pts eq V`pts then
            for j in [1..#V`pts] do
                if &and[ v eq 0 : v in V`vals[j] ] then V`vals[j] := U`vals[j]; end if;
            end for;
        else
            nc := V`ncirc;
            V`pts := V`pts cat U`pts; V`owner := V`owner cat [ o + nc : o in U`owner ];
            V`isctr := V`isctr cat U`isctr; V`vals := V`vals cat U`vals; V`ncirc := nc + U`ncirc;
        end if;
    end for;
    return V;
end intrinsic;

intrinsic WriteNelsonValues(V::Rec, filename::MonStgElt)
{Write a value record (NelsonCircleValues) as a Magma file that ReadNelsonValues reads back.}
    F := Open(filename, "w");
    fprintf F, "// Nelson weight-2 values for X_0^%o(%o)/%o, cut %o, prec %o, rho %o\n",
            V`D, V`N, V`W, V`cut, V`prec, V`rho;
    fprintf F, "DD := %o; NN := %o; WW := %o; GG := %o; PRECV := %o; CUTV := %o;\n",
            V`D, V`N, V`W, V`g, V`prec, V`cut;
    fprintf F, "RHO := %o; NCIRC := %o; NPC := %o;\n", V`rho, V`ncirc, V`npc;
    fprintf F, "HALLS := %o;\nALSIG := %o;\nOBAS := %o;\n", V`halls, V`alsig, V`obas;
    fprintf F, "PTS := [\n";
    for j in [1..#V`pts] do
        fprintf F, "  [%o, %o]%o\n", Real(V`pts[j]), Imaginary(V`pts[j]), j lt #V`pts select "," else "";
    end for;
    fprintf F, "];\nOWNER := %o;\nISCTR := %o;\n", V`owner, V`isctr;
    fprintf F, "VALS := [\n";
    for j in [1..#V`vals] do
        fprintf F, "  [%o]%o\n",
            &cat[ Sprintf("[%o, %o]%o", Real(V`vals[j][i]), Imaginary(V`vals[j][i]), i lt V`g select ", " else "")
                  : i in [1..V`g] ],
            j lt #V`vals select "," else "";
    end for;
    fprintf F, "];\n";
    delete F;
end intrinsic;

intrinsic ReadNelsonValues(filename::MonStgElt) -> Rec
{Read a value record written by WriteNelsonValues.}
    txt := Read(filename);
    hasobas := Position(txt, "OBAS") gt 0;
    dmp := eval (txt cat "\nreturn <DD, NN, WW, GG, PRECV, CUTV, RHO, NCIRC, NPC, HALLS, ALSIG, PTS, OWNER, ISCTR, VALS" cat
                 (hasobas select ", OBAS" else "") cat ">;");
    CC := ComplexField(dmp[5]); Ic := CC.1;
    pts := [ CC | p[1] + Ic*p[2] : p in dmp[12] ];
    vals := [ [ CC | v[1] + Ic*v[2] : v in row ] : row in dmp[15] ];
    return rec< NelVals | D := dmp[1], N := dmp[2], W := dmp[3], g := dmp[4], prec := dmp[5], cut := dmp[6],
                          rho := RealField(dmp[5]) ! dmp[7], ncirc := dmp[8], npc := dmp[9], halls := dmp[10],
                          alsig := dmp[11], pts := pts, owner := dmp[13], isctr := dmp[14], vals := vals,
                          obas := hasobas select dmp[16] else [] >;
end intrinsic;

/////////////////////////////////////////////////////////////////////////////////////////////////
// The model.  An Atkin-Lehner involution u with X/<W,u> of genus 1 fixes the coordinate on the
// canonical conic (u : x -> -x), which makes the octic even and its quartic in t = x^2 carry the
// three invariants c0c2/c1^2, c1c3/c2^2, c2c4/c3^2 -- of small enough height to recognise, unlike
// the coordinate-free invariants of the octic.  Both the evenness and the u-equivariance of the
// conic are read off the fit, not assumed.  See paper/nelson-method.tex, "The model".

// the classical certainty criterion: a rational p/q within eps of a real number is accepted when
// q^2 eps < 1/2, so that no other fraction of denominator at most q lies as close
function RecogniseRational(z, acc)
    r := Real(z);
    scale := Abs(r) gt 1 select Abs(r) else Parent(r) ! 1;
    eps := Max(10^(-acc) * scale, 20 * Abs(Imaginary(z)));
    for e in [2..acc] do
        q := BestApproximation(r, 10^e);
        if Abs(r - q) lt eps and Denominator(q)^2 * eps lt 1/2 then return true, q; end if;
    end for;
    return false, 0;
end function;

// real SVD of a complex matrix: singular values (ascending, once each) and a kernel vector
function ComplexKernel(My)
    n2 := Nrows(My); nm := Ncols(My);
    CC := BaseRing(My); RR := RealField(Precision(CC)); Ic := CC.1;
    MR := Matrix(RR, 2*n2, 2*nm, [ RR ! 0 : i in [1..4*n2*nm] ]);
    for i in [1..n2] do for j in [1..nm] do
        a := Real(My[i][j]); b := Imaginary(My[i][j]);
        MR[i][j] := a; MR[i][nm+j] := -b; MR[n2+i][j] := b; MR[n2+i][nm+j] := a;
    end for; end for;
    Sg, U, V := SingularValueDecomposition(MR);
    rr := Min(Nrows(MR), Ncols(MR));
    sv := Reverse([ Abs(Sg[i][i]) : i in [1..rr] ]);
    sv := [ sv[i] : i in [1..#sv by 2] ];
    // which factor holds the right singular vectors, and as rows or columns, is settled by
    // trying both readings and keeping the one that annihilates the matrix
    kcol := Ncols(V);
    function KernRes(vec)
        s := RR ! 0;
        for i in [1..n2] do s +:= Abs(&+[ My[i][c]*vec[c] : c in [1..nm] ])^2; end for;
        return Sqrt(s)/Sqrt(&+[ Abs(c)^2 : c in vec ]);
    end function;
    vcol := [ CC | V[j][kcol] + Ic*V[nm+j][kcol] : j in [1..nm] ];
    vrow := [ CC | V[kcol][j] + Ic*V[kcol][nm+j] : j in [1..nm] ];
    rc := KernRes(vcol); rw := KernRes(vrow);
    return sv, (rc lt rw select vcol else vrow), Minimum(rc, rw);
end function;

// The geometry of the model step -- conic, coordinate, Cauchy derivative, quartic -- which
// NelsonMarkedInvariants exposes and NelsonHyperellipticModel continues by recognising the
// three invariants.  Returns them as complex numbers together with the involution and the
// diagnostics.
function MarkedInvariantsOf(V, Involution)
    g := V`g;
    error if not (g eq 3), "genus 3 only";
    prec := V`prec;
    CC := ComplexField(prec); RR := RealField(prec); Ic := CC.1; R6 := RealField(8);
    np := V`npc; nc := V`ncirc; rho := RR ! V`rho;
    pts := V`pts; HALLS := V`halls; ALSIG := V`alsig; ISCTR := V`isctr;
    error if not (#pts eq nc*(np+1)), "point list does not tile the circles";
    error if not (forall{ j : j in [1..#pts] | exists{ v : v in V`vals[j] | v ne 0 } }), "some points carry no value (incomplete chunks?)";
    Gv := [ [ CC | V`vals[j][i] / Imaginary(pts[j]) : i in [1..g] ] : j in [1..#pts] ];
    diag := rec< recformat< u, conicsv, conicres, equivariance, normalform, cauchy, centremean, xscale,
                            worstprediction, oddcoefficients, lsresidual, A, B, C, Arat, Brat, Crat,
                            recognised, kappa, mu, quartic, quotient, counts > | >;

    // the involution
    uix := 0;
    if Involution ne 0 then
        uix := Index(HALLS, Involution);
        error if not (uix ne 0), "Involution is not a Hall divisor of D*N";
        error if not (#[ r : r in ALSIG | r[uix] eq 1 ] eq 1), "X/<W,u> does not have genus 1 for this u";
    else
        for i in [1..#HALLS] do
            if HALLS[i] eq 1 then continue; end if;
            if #[ r : r in ALSIG | r[i] eq 1 ] eq 1 then uix := i; break; end if;
        end for;
        error if not (uix ne 0), "no Atkin-Lehner involution with a genus-1 quotient";
    end if;
    um := HALLS[uix]; diag`u := um;
    pind := [ j : j in [1..g] | ALSIG[j][uix] eq 1 ][1];
    aind := [ j : j in [1..g] | ALSIG[j][uix] eq -1 ];

    // the conic
    idx := [ j : j in [1..#pts] | ISCTR[j] eq 1 ] cat [ j : j in [1..#pts] | ISCTR[j] eq 0 and (j mod np) eq 3 ];
    mons := [ [a,b] : b in [a..g], a in [1..g] ];
    Mx := Matrix(CC, #idx, #mons, &cat[ [ Gv[j][m[1]]*Gv[j][m[2]] : m in mons ] : j in idx ]);
    for i in [1..#idx] do
        nm := Sqrt(&+[ Abs(Mx[i][c])^2 : c in [1..#mons] ]);
        for c in [1..#mons] do Mx[i][c] := Mx[i][c]/nm; end for;
    end for;
    cn := [ Sqrt(&+[ Abs(Mx[i][c])^2 : i in [1..#idx] ]) : c in [1..#mons] ];
    for c in [1..#mons] do for i in [1..#idx] do Mx[i][c] := Mx[i][c]/cn[c]; end for; end for;
    sv, kvec, kres := ComplexKernel(Mx);
    diag`conicsv := sv; diag`conicres := kres;
    vprintf NelsonShimizu, 1: "conic singular values %o\n", [ R6 ! x : x in sv ];
    error if not (sv[2]/sv[1] ge 10^4), "no clean conic through the canonical image";
    error if not (kres le 10^(-8)), "neither reading of the SVD gives a kernel vector";
    qc := [ CC | kvec[j]/cn[j] : j in [1..#mons] ];
    Q := ZeroMatrix(CC, 3, 3);
    for c in [1..#mons] do
        a := mons[c][1]; b := mons[c][2];
        if a eq b then Q[a][a] := qc[c]; else Q[a][b] := qc[c]/2; Q[b][a] := qc[c]/2; end if;
    end for;
    scl := Maximum([ Abs(Q[a][b]) : a, b in [1..3] ]);
    diag`equivariance := Maximum(Abs(Q[pind][aind[1]]), Abs(Q[pind][aind[2]]))/scl;
    vprintf NelsonShimizu, 1: "u = w_%o; cross terms of the conic between the invariant and anti-invariant differentials %o\n",
        um, R6 ! diag`equivariance;

    // the coordinate: Q = qpp ep^2 + (binary form in ea, eb), ep^2 = (-qaa/qpp)(ea - r1 eb)(ea - r2 eb)
    a1 := aind[1]; a2 := aind[2];
    qpp := Q[pind][pind]; qaa := Q[a1][a1]; qab := Q[a1][a2]; qbb := Q[a2][a2];
    disc := qab^2 - qaa*qbb;
    r1 := (-qab + Sqrt(disc))/qaa; r2 := (-qab - Sqrt(disc))/qaa;
    kap := -qaa/qpp;
    function Wof(j)
        ep := Gv[j][pind]; ea := Gv[j][a1]; eb := Gv[j][a2];
        return ea - r1*eb, ep, kap*(ea - r2*eb);
    end function;
    nf := RR ! 0;
    for j in [1..#pts] do
        w1, w2, w3 := Wof(j);
        nf := Max(nf, Abs(w2^2 - w1*w3)/(Abs(w2)^2 + Abs(w1*w3)));
    end for;
    diag`normalform := nf;

    // x and dx/dz at the centres
    xs := [ CC | ]; ys := [ CC | ]; dchk := RR ! 0; cchk := RR ! 0;
    for c in [1..nc] do
        base := (c-1)*(np+1) + 1;
        xc := [ CC | ];
        for t in [1..np] do
            w1, w2, w3 := Wof(base+t);
            Append(~xc, w2/w1);
        end for;
        w1c, w2c, w3c := Wof(base);
        xcen := w2c/w1c;
        om := Exp(2*Pi(CC)*Ic/np);
        d1 := (&+[ xc[t+1]*om^(-t) : t in [0..np-1] ])/(np*rho);
        d0 := (&+[ xc[t+1] : t in [0..np-1] ])/np;
        d1h := (&+[ xc[2*t+1]*om^(-2*t) : t in [0..(np div 2)-1] ])/((np div 2)*rho);
        dchk := Max(dchk, Abs(d1-d1h)/Abs(d1));
        cchk := Max(cchk, Abs(d0-xcen)/Abs(xcen));
        Append(~xs, xcen); Append(~ys, d1/w1c);
    end for;
    diag`cauchy := dchk; diag`centremean := cchk;
    vprintf NelsonShimizu, 1: "derivative by %o against %o circle points: %o; |x| at the centres %o\n",
        np, np div 2, R6 ! dchk, [ R6 ! Abs(x) : x in xs ];
    // x -> lambda x is the one freedom u leaves; rescaling to modulus about 1 conditions the
    // Vandermonde system and leaves A, B, C unchanged
    sigx := Exp((&+[ Log(Abs(xx)) : xx in xs ])/nc);
    diag`xscale := sigx;
    xs := [ xx/sigx : xx in xs ]; ys := [ yy/sigx : yy in ys ];

    // the quartic in t = x^2: five unknowns, the other centres predictions
    ts := [ xx^2 : xx in xs ];
    nun := 5;
    error if not (nc ge nun + 1), "at least six circles are needed to leave a prediction";
    Amx := Matrix(CC, nun, nun, &cat[ [ ts[j]^(i-1) : i in [1..nun] ] : j in [1..nun] ]);
    bv := Matrix(CC, nun, 1, [ ys[j]^2 : j in [1..nun] ]);
    cf := Amx^(-1)*bv;
    coef := [ CC | cf[i][1] : i in [1..nun] ];
    maxpr := RR ! 0;
    for j in [nun+1..nc] do
        pv := &+[ coef[i]*ts[j]^(i-1) : i in [1..nun] ];
        maxpr := Maximum(maxpr, Abs(pv - ys[j]^2)/Abs(ys[j]^2));
    end for;
    diag`worstprediction := maxpr;
    if nc ge 10 then
        Ao := Matrix(CC, 9, 9, &cat[ [ xs[j]^(i-1) : i in [1..9] ] : j in [1..9] ]);
        bo := Matrix(CC, 9, 1, [ ys[j]^2 : j in [1..9] ]);
        co := Ao^(-1)*bo;
        oc := [ CC | co[i][1] : i in [1..9] ];
        sco := Maximum([ Abs(c) : c in oc ]);
        diag`oddcoefficients := [ RR | Abs(oc[i])/sco : i in [2,4,6,8] ];
    end if;
    // the values from every centre, each equation divided by y^2 so that no centre dominates
    Aw := Matrix(CC, nc, nun, &cat[ [ ts[j]^(i-1)/ys[j]^2 : i in [1..nun] ] : j in [1..nc] ]);
    bw := Matrix(CC, nc, 1, [ CC | 1 : j in [1..nc] ]);
    Ah := Matrix(CC, nun, nc, &cat[ [ ComplexConjugate(Aw[j][i]) : j in [1..nc] ] : i in [1..nun] ]);
    cls := (Ah*Aw)^(-1)*(Ah*bw);
    gc := [ CC | cls[i][1] : i in [1..nun] ];
    diag`lsresidual := Maximum([ RR | Abs(&+[ gc[i]*ts[j]^(i-1) : i in [1..nun] ]/ys[j]^2 - 1) : j in [1..nc] ]);
    vprintf NelsonShimizu, 1: "quartic in x^2: worst prediction %o, least-squares residual %o\n",
        R6 ! maxpr, R6 ! diag`lsresidual;

    A := gc[1]*gc[3]/gc[2]^2; Bq := gc[2]*gc[4]/gc[3]^2; Cq := gc[3]*gc[5]/gc[4]^2;
    diag`A := A; diag`B := Bq; diag`C := Cq;
    return A, Bq, Cq, um, diag;
end function;

intrinsic NelsonMarkedInvariants(V::Rec : Involution := 0) -> FldComElt, FldComElt, FldComElt, RngIntElt, Rec
{The three marked invariants A, B, C of the quartic of a genus-3 value record, as complex numbers,
 with the Atkin-Lehner element u that fixes the coordinate and the diagnostics of the fit (conic
 singular values, u-equivariance, derivative convergence, worst prediction, least-squares
 residual).  This is NelsonHyperellipticModel without the step that recognises them as rational
 numbers, which is what a comparison against a model known by other means needs: the invariants
 of such a model can be compared with these whether or not the certainty criterion accepts them.
 Numerical.}
    A, Bq, Cq, um, diag := MarkedInvariantsOf(V, Involution);
    return A, Bq, Cq, um, diag;
end intrinsic;

intrinsic NelsonHyperellipticModel(V::Rec : Involution := 0, Digits := 0, PrimeCap := 150) -> CrvHyp, Rec
{From a value record of a genus-3 quotient (NelsonCircleValues, or ReadNelsonValues) to a
 hyperelliptic model over Q.  Involution names the Atkin-Lehner element u pinning the coordinate
 (default: the first Hall divisor outside W under which exactly one column is invariant, i.e.
 X/<W,u> has genus 1).  Steps: the canonical conic from the centres and one point of each circle;
 the coordinate x with u : x -> -x; x and dx/dz at each centre by Cauchy's formula; the quartic
 g in t = x^2 by least squares, with the centres beyond the first five as predictions; the
 invariants A, B, C recognised by the certainty criterion with the noise taken from their
 imaginary parts (Digits caps the digits claimed; default a third of the precision); then
 NelsonTwistConstants.  When B is rational but A and C are not, the fixed points of u on the
 conic are conjugate over a quadratic field K; K is read off the fields of definition of the
 fixed points of u on X, A is recognised in K, and NelsonQuadraticTwist gives the curve over Q
 (or, when the conic has no rational point, the K-model with its descent datum, returned as a
 curve over K; kappa and mu then hold r and e).  Fails with a message at the first step that
 does not resolve.  The record returned carries every diagnostic (conic singular values,
 u-equivariance, derivative convergence, worst prediction, the complex invariants, kappa, mu, the
 quartic).  Numerical: the model is recognised, not derived, and the point counts confirm only
 the isogeny class.}
    A, Bq, Cq, um, diag := MarkedInvariantsOf(V, Involution);
    g := V`g; prec := V`prec;
    CC := ComplexField(prec); RR := RealField(prec); R6 := RealField(8);
    acc := Digits eq 0 select Maximum(8, prec div 3) else Digits;
    okA, Ar := RecogniseRational(A, acc);
    okB, Br := RecogniseRational(Bq, acc);
    okC, Cr := RecogniseRational(Cq, acc);
    diag`recognised := [okA, okB, okC];
    vprintf NelsonShimizu, 1: "invariants (noise %o): A = %o, B = %o, C = %o\n",
        [ R6 ! Abs(Imaginary(z)) : z in [A, Bq, Cq] ],
        okA select Sprint(Ar) else "not recognised", okB select Sprint(Br) else "not recognised",
        okC select Sprint(Cr) else "not recognised";
    if okB and not (okA and okC) then
        // B rational with A, C irrational is the signature of fixed points of u on the conic that
        // are conjugate over a quadratic field K: K comes from the fields of definition of the
        // fixed points of u on X, A is recognised in K, and the curve from NelsonQuadraticTwist.
        dK, info := NelsonBranchField(V`D, V`N, Set(V`W), um);
        require dK ne 0 : "B is rational but A, C are not, and the fixed points of u on X show no common quadratic field";
        dsf := Squarefree(dK);
        noise := Maximum([ RR | Abs(Imaginary(Bq)), Abs(Real(A) - Real(Cq)), Abs(Imaginary(A) + Imaginary(Cq))/2,
                           dsf gt 0 select Abs(Imaginary(A)) else 0, 10^(-acc) ]);
        okK, AK, bound := NelsonRecogniseInQuadraticField(A, Cq, dsf, RR ! (2*noise));
        vprintf NelsonShimizu, 1: "fixed points of u conjugate over Q(sqrt %o); noise %o, so denominators are certified to %o\n",
            dsf, RealField(6) ! noise, bound;
        require okK : Sprintf("the fixed points of u are conjugate over Q(sqrt %o) and B = %o, but A is not recognised in that field at accuracy %o (denominators certified up to %o)", dsf, Br, RealField(4) ! noise, bound);
        diag`Arat := AK; diag`Brat := Br; diag`Crat := Parent(AK) ! 0;
        diag`recognised := [true, true, true];
        vprintf NelsonShimizu, 1: "A = %o in Q(sqrt %o)\n", AK, dsf;
        e, r, R := NelsonQuadraticTwist(V`D, V`N, Set(V`W), um, AK, Br : PrimeCap := PrimeCap);
        diag`kappa := r; diag`mu := e; diag`quotient := R`quotient; diag`counts := R`counts;
        diag`quartic := R`PE;
        if R`split then return R`model, diag; end if;
        return HyperellipticCurve(R`F), diag;
    end if;
    require okA and okB and okC :
        Sprintf("invariants not recognised at the certainty criterion: A %o, B %o, C %o",
                okA select Sprint(Ar) else "?", okB select Sprint(Br) else "?", okC select Sprint(Cr) else "?");
    diag`Arat := Ar; diag`Brat := Br; diag`Crat := Cr;
    kappa, mu, Cfin, Ecur, gnorm, ncounts := NelsonTwistConstants(V`D, V`N, Set(V`W), um, [Ar, Br, Cr] : PrimeCap := PrimeCap);
    diag`kappa := kappa; diag`mu := mu; diag`quartic := gnorm; diag`quotient := Ecur; diag`counts := ncounts;
    return Cfin, diag;
end intrinsic;

intrinsic NelsonTwistConstants(D::RngIntElt, N::RngIntElt, W::SetEnum, um::RngIntElt,
                               ABC::SeqEnum : PrimeCap := 150)
          -> RngIntElt, RngIntElt, CrvHyp, CrvHyp, RngUPolElt, RngIntElt
{Given the three marked invariants [A, B, C] of the quartic g, so that the curve is y^2 = kappa g0(mu x^2)
 with g0(t) = 1 + t + A t^2 + A^2 B t^3 + A^3 B^2 C t^4 and kappa, mu rational numbers that matter
 only modulo squares (a complex-analytic computation cannot see them), the two numbers solved
 rather than searched: each odd prime p of good reduction compares the point count of
 Y^2 = g0(t) over F_p with the trace formula for X/<W,u> and gives the Legendre symbol (kappa|p),
 one linear equation over F_2 in the exponents of kappa over the primes of 2DN and the bad primes
 of g0; then mu likewise from the curve itself.  Both systems must have a unique solution.
 Returns kappa, mu, the reduced minimal model of the curve, the genus-1 quotient Y^2 = kappa g0,
 g0, and the number of point counts over F_p and F_p^2 (p < PrimeCap) that the model matched,
 which fixes the isogeny class of the Jacobian and nothing more.  The invariants are numerical
 input recognised as rationals; the counts are exact.}
    require #ABC eq 3 : "three invariants expected";
    A := Rationals() ! ABC[1]; B := Rationals() ! ABC[2]; C := Rationals() ! ABC[3];
    PQ<x> := PolynomialRing(Rationals());
    gnorm := 1 + x + A*x^2 + A^2*B*x^3 + A^3*B^2*C*x^4;
    require Degree(gnorm) eq 4 and Discriminant(gnorm) ne 0 : "the quartic is degenerate";
    DN := D*N;
    WU := { ALMul(a, um) : a in W } join W;
    gE := GenusShimuraCurveQuotient(D, N, WU);
    gX := GenusShimuraCurveQuotient(D, N, W);
    require gX eq 3 and gE eq 1 : "expected genus 3 over a genus-1 Atkin-Lehner quotient";
    XQ := CreateShimuraQuot(D, N, W);  XQ`g := gX;
    EQ := CreateShimuraQuot(D, N, WU); EQ`g := gE;
    good := [ p : p in PrimesUpTo(PrimeCap) | p ne 2 and DN mod p ne 0 ];
    tX := AssociativeArray(); tE := AssociativeArray();
    for p in good do
        tX[p] := [ ComputePointsViaTrace(XQ, p, d) : d in [1,2] ];
        tE[p] := [ ComputePointsViaTrace(EQ, p, d) : d in [1,2] ];
    end for;
    // support of the two unknowns: the primes of 2DN and the bad primes of the quartic
    cs := Coefficients(gnorm);
    dd := LCM([ Denominator(c) : c in cs ]);
    S := Set(PrimeDivisors(2*DN)) join Set(PrimeDivisors(Integers() ! (dd^4 * Discriminant(gnorm))));
    for c in cs do
        if c ne 0 then S join:= Set(PrimeDivisors(Numerator(c))); end if;
        S join:= Set(PrimeDivisors(Denominator(c)));
    end for;
    supp := Sort(SetToSequence(S));
    // a prime dividing a denominator is reported as bad rather than counted on a model it does
    // not belong to; scaling by the square of the common denominator is an isomorphism over Q
    function CountQ(poly, p, d)
        den := LCM([ Denominator(c) : c in Coefficients(poly) ] cat [1]);
        if den mod p eq 0 then return -1; end if;
        pint := PQ ! (den^2 * poly);
        gg := ChangeRing(pint, GF(p^d));
        if Degree(gg) lt Degree(pint) - 1 or Degree(gg) lt 3 or Discriminant(gg) eq 0 then return -1; end if;
        return #Points(HyperellipticCurve(gg));
    end function;
    function SolveSymbols(rows, rhs, name)
        error if #rows lt #supp + 4, "too few primes for " cat name cat "; raise PrimeCap";
        M := Matrix(GF(2), #rows, #supp+1, &cat rows);
        ok, v := IsConsistent(Transpose(M), Vector(GF(2), rhs));
        error if not ok, "the Legendre symbols for " cat name cat " are inconsistent: the quartic is not a twist of the quotient";
        error if Rank(M) ne #supp+1, name cat " is not determined uniquely; raise PrimeCap";
        sgn := (v[#supp+1] eq GF(2)!1) select -1 else 1;
        return sgn * (&*([Integers() | 1] cat [ supp[i] : i in [1..#supp] | v[i] eq GF(2)!1 ]));
    end function;
    function SymbolRow(p)
        return [ GF(2) | (LegendreSymbol(q, p) eq 1 select 0 else 1) : q in supp ]
               cat [ GF(2) | (p mod 4 eq 1 select 0 else 1) ];
    end function;
    rows := []; rhs := [];
    for p in good do
        n1 := CountQ(gnorm, p, 1);
        if n1 lt 0 then continue; end if;
        t1 := tE[p][1];
        if n1 eq t1 and 2*(p+1) - n1 eq t1 then continue; end if;
        if n1 eq t1 then b := 0; elif 2*(p+1) - n1 eq t1 then b := 1; else continue; end if;
        Append(~rows, SymbolRow(p)); Append(~rhs, GF(2) ! b);
    end for;
    kappa := SolveSymbols(rows, rhs, "kappa");
    Ecur := HyperellipticCurve(kappa*gnorm);
    for p in good do for d in [1,2] do
        c := CountQ(kappa*gnorm, p, d);
        require c lt 0 or c eq tE[p][d] : Sprintf("the quotient's count disagrees at p = %o, degree %o", p, d);
    end for; end for;
    rows := []; rhs := [];
    for p in good do
        f0 := kappa*Evaluate(gnorm, x^2);
        nsq := 1; while LegendreSymbol(nsq, p) eq 1 do nsq +:= 1; end while;
        f1 := kappa*Evaluate(gnorm, nsq*x^2);
        n0 := CountQ(f0, p, 1); n1 := CountQ(f1, p, 1);
        if n0 lt 0 or n1 lt 0 or n0 eq n1 then continue; end if;
        t1 := tX[p][1];
        if n0 eq t1 then b := 0; elif n1 eq t1 then b := 1; else continue; end if;
        Append(~rows, SymbolRow(p)); Append(~rhs, GF(2) ! b);
    end for;
    mu := SolveSymbols(rows, rhs, "mu");
    foct := kappa*Evaluate(gnorm, mu*x^2);
    dd := LCM([ Denominator(c) : c in Coefficients(foct) ]);
    Cfin := HyperellipticCurve(dd^2*foct);
    try Cfin := ReducedMinimalWeierstrassModel(Cfin); catch e
        try Cfin := ReducedModel(SimplifiedModel(Cfin)); catch e2 ; end try;
    end try;
    require Genus(Cfin) eq 3 : "genus mismatch";
    fq, hq := HyperellipticPolynomials(Cfin);
    gg := fq + hq^2/4;
    nchk := 0;
    for p in good do for d in [1,2] do
        c := CountQ(gg, p, d);
        if c lt 0 then continue; end if;
        require c eq tX[p][d] : Sprintf("point count disagrees at p = %o, degree %o: %o against %o", p, d, c, tX[p][d]);
        nchk +:= 1;
    end for; end for;
    vprintf NelsonShimizu, 1: "kappa = %o, mu = %o; %o point counts agree; model %o\n", kappa, mu, nchk, Cfin;
    return kappa, mu, Cfin, Ecur, gnorm, nchk;
end intrinsic;

intrinsic NelsonHyperellipticModel(D::RngIntElt, N::RngIntElt, W::SetEnum :
                                   Precision := 80, Cutoff := 46, Circles := 14, CirclePoints := 16,
                                   Radius := 3/100, Centres := [], Involution := 0, Digits := 0,
                                   PrimeCap := 150, OrderBasis := [], ChooseBasePoint := true,
                                   AllowCoarsePrecision := false) -> CrvHyp, Rec, Rec
{The whole route from (D, N, W), genus 3 with an Atkin-Lehner involution u such that X/<W,u> has
 genus 1, to a hyperelliptic model over Q: Circles centres (NelsonSamplePoints, or Centres when
 given) with CirclePoints points on a circle of radius Radius around each, evaluated by Theorem
 3.1 at Cutoff and Precision, then NelsonHyperellipticModel on the values.  Returns the model, the
 diagnostics record and the value record (for WriteNelsonValues).  The cost is Circles *
 (CirclePoints + 1) evaluations, each linear in DN; at cutoff 46 the precision should be 80 digits
 or more.  Numerical.}
    B := NelsonBasis(D, N, W : Precision := Precision, Cutoff := Cutoff, OrderBasis := OrderBasis,
                              AllowCoarsePrecision := AllowCoarsePrecision);
    require B`g eq 3 : "genus 3 only";
    if ChooseBasePoint then NelsonChooseBasePoint(~B); end if;
    ctr := #Centres eq 0 select NelsonSamplePoints(B, Circles) else [ B`CC ! c : c in Centres ];
    V := NelsonCircleValues(B, ctr, CirclePoints, B`RR ! Radius);
    C, diag := NelsonHyperellipticModel(V : Involution := Involution, Digits := Digits, PrimeCap := PrimeCap);
    return C, diag, V;
end intrinsic;

intrinsic NelsonCoordinateSpread(V::Rec : Involution := 0) -> SeqEnum
{The moduli |x| of the coordinate x = w2/w1 (as in NelsonHyperellipticModel) at the centres of a
 value record, relative to their geometric mean.  The quartic in t = x^2 is determined at its
 large-t end by the centres with large |x| and at its small-t end by the others; when the moduli
 span less than a factor of about five the fit is ill-conditioned and the invariants come out to
 far fewer digits than the values carry, so centres should be added where |x| is large (a coarse
 evaluation at a small cutoff is enough to find them, since the ratio is read off the conic).  Numerical.}
    g := V`g;
    prec := V`prec; CC := ComplexField(prec); RR := RealField(prec);
    np := V`npc; pts := V`pts; HALLS := V`halls; ALSIG := V`alsig; ISCTR := V`isctr;
    Gv := [ [ CC | V`vals[j][i] / Imaginary(pts[j]) : i in [1..g] ] : j in [1..#pts] ];
    uix := 0;
    for i in [1..#HALLS] do
        if HALLS[i] eq 1 then continue; end if;
        if Involution ne 0 and HALLS[i] ne Involution then continue; end if;
        if #[ r : r in ALSIG | r[i] eq 1 ] eq 1 then uix := i; break; end if;
    end for;
    require uix ne 0 : "no suitable involution";
    pind := [ j : j in [1..g] | ALSIG[j][uix] eq 1 ][1];
    aind := [ j : j in [1..g] | ALSIG[j][uix] eq -1 ];
    idx := [ j : j in [1..#pts] | ISCTR[j] eq 1 ] cat [ j : j in [1..#pts] | ISCTR[j] eq 0 and (j mod np) eq 3 ];
    mons := [ [a,b] : b in [a..g], a in [1..g] ];
    Mx := Matrix(CC, #idx, #mons, &cat[ [ Gv[j][m[1]]*Gv[j][m[2]] : m in mons ] : j in idx ]);
    for i in [1..#idx] do
        nm := Sqrt(&+[ Abs(Mx[i][c])^2 : c in [1..#mons] ]);
        for c in [1..#mons] do Mx[i][c] := Mx[i][c]/nm; end for;
    end for;
    cn := [ Sqrt(&+[ Abs(Mx[i][c])^2 : i in [1..#idx] ]) : c in [1..#mons] ];
    for c in [1..#mons] do for i in [1..#idx] do Mx[i][c] := Mx[i][c]/cn[c]; end for; end for;
    sv, kvec, kres := ComplexKernel(Mx);
    require kres le 10^(-5) : "no conic";
    qc := [ CC | kvec[j]/cn[j] : j in [1..#mons] ];
    Q := ZeroMatrix(CC, 3, 3);
    for c in [1..#mons] do
        a := mons[c][1]; b := mons[c][2];
        if a eq b then Q[a][a] := qc[c]; else Q[a][b] := qc[c]/2; Q[b][a] := qc[c]/2; end if;
    end for;
    a1 := aind[1]; a2 := aind[2];
    qaa := Q[a1][a1]; qab := Q[a1][a2]; qbb := Q[a2][a2];
    r1 := (-qab + Sqrt(qab^2 - qaa*qbb))/qaa;
    xs := [ RR | Abs(Gv[j][pind]/(Gv[j][a1] - r1*Gv[j][a2])) : j in [1..#pts] | ISCTR[j] eq 1 ];
    sigx := Exp((&+[ Log(x) : x in xs ])/#xs);
    return [ x/sigx : x in xs ];
end intrinsic;

/////////////////////////////////////////////////////////////////////////////////////////////////
// Descent from a quadratic field.  When the two fixed points of u on the canonical conic are
// conjugate over K = Q(sqrt d), the quartic lives over K (B rational, C = sigma(A)), and the
// Galois-stability of the branch divisor is a norm equation in K, solvable exactly when the
// K-model is the base change of a curve over Q with this marking (Hilbert 90).  Two twists remain,
// the twist by u and the quadratic twist, both read off point counts; when the twist by u is not a
// norm from K the conic X/iota has no rational point and no model y^2 = f(s) over Q exists.
// See paper/nelson-method.tex, "When the coordinate is not rational".

// a Legendre-symbol row for the F_2 system in the exponents of a twist constant over supp
function TwistSymbolRow(supp, p)
    return [ GF(2) | (LegendreSymbol(q, p) eq 1 select 0 else 1) : q in supp ]
           cat [ GF(2) | (p mod 4 eq 1 select 0 else 1) ];
end function;

function TwistSolve(supp, rows, rhs, name)
    error if #rows lt #supp + 4, "too few primes for " cat name cat "; raise PrimeCap";
    M := Matrix(GF(2), #rows, #supp+1, &cat rows);
    ok, v := IsConsistent(Transpose(M), Vector(GF(2), rhs));
    error if not ok, "the Legendre symbols for " cat name cat " are inconsistent: the curve is not a twist of the target";
    error if Rank(M) ne #supp+1, name cat " is not determined uniquely; raise PrimeCap";
    sgn := (v[#supp+1] eq GF(2)!1) select -1 else 1;
    return sgn * (&*([Integers() | 1] cat [ supp[i] : i in [1..#supp] | v[i] eq GF(2)!1 ]));
end function;

// point count of y^2 = poly over F_{p^d}, or -1 when the model is bad there
function TwistCountQ(poly, p, d)
    PQ := Parent(poly);
    den := LCM([ Denominator(c) : c in Coefficients(poly) ] cat [1]);
    if den mod p eq 0 then return -1; end if;
    pint := PQ ! (den^2 * poly);
    gg := ChangeRing(pint, GF(p^d));
    if Degree(gg) lt Degree(pint) - 1 or Degree(gg) lt 3 or Discriminant(gg) eq 0 then return -1; end if;
    return #Points(HyperellipticCurve(gg));
end function;

function PolySupport(poly, extra)
    cs := Coefficients(poly);
    dd := LCM([ Denominator(c) : c in cs ]);
    S := Set(PrimeDivisors(extra)) join Set(PrimeDivisors(Integers() ! (dd^(2*Degree(poly)) * Discriminant(poly))));
    for c in cs do
        if c ne 0 then S join:= Set(PrimeDivisors(Numerator(c))); end if;
        S join:= Set(PrimeDivisors(Denominator(c)));
    end for;
    return Sort(SetToSequence(S));
end function;

intrinsic NelsonQuadraticDescent(A::FldNumElt, B::FldRatElt : Twist := 1) -> RngUPolElt, RngUPolElt, FldNumElt, RngUPolElt
{Given the marked invariants of the quartic over a quadratic field K = Q(sqrt d), A in K and B
 rational (C = sigma(A)), of a genus-3 curve whose marking involution has conjugate fixed points
 on the canonical conic: solves the norm equation N(mu) = 1/(B N(A)) in K (fails, with a message,
 when 1/(B N(A)) is not a norm, in which case the K-model is not the base change of a curve over
 Q with this marking), forms F(x) = G0(mu x^2) over K, checks the Galois condition
 sigma(F)(x) = c x^8 F(1/x), and substitutes x = (1 + w s)/(1 - w s).  Returns P(s) in Q[s]
 (the branch octic over Q, defined up to a rational constant: the curve is y^2 = r P(s)), the
 quartic P_E(v) in Q[v] of the quotient by u in the same way (t = x^2 = (1 + w v)/(1 - w v);
 the quotient is y^2 = r_E P_E(v)), mu, and F.  The solution mu of the norm equation is defined
 up to norm-one elements sigma(gamma)/gamma, and replacing mu by mu sigma(gamma)/gamma gives the
 twist of the curve by the involution u over Q(sqrt N(gamma)): the forms of the marked curve over
 Q with a split conic are the u-twists by the rational numbers e that are norms from K, together
 with the quadratic twists by r.  Twist := e (a norm from K) selects that form; which e is the
 curve's is a matter of point counts (NelsonQuadraticTwist).  The quotient by u is the same for
 every e.  Exact.}
    K := Parent(A);
    require Degree(K) eq 2 : "A must lie in a quadratic field";
    w := K.1;
    d := Rationals() ! (w^2);
    require d in Rationals() : "the generator of K must be a square root of a rational number";
    sig := hom< K -> K | -w >;
    C := sig(A);
    PK<x> := PolynomialRing(K);
    G0 := 1 + x + A*x^2 + A^2*B*x^3 + A^3*B^2*C*x^4;
    require Degree(G0) eq 4 and Discriminant(G0) ne 0 : "the quartic is degenerate";
    target := 1/(B*Norm(A));
    ok, sols := NormEquation(K, target);
    require ok : Sprintf("1/(B N(A)) = %o is not a norm from %o: the quartic with this marking does not descend", target, K);
    mu := sols[1];
    e := Rationals() ! Twist;
    if e ne 1 then
        oke, gs := NormEquation(K, e);
        require oke : Sprintf("the twist %o is not a norm from %o", e, K);
        mu := mu * sig(gs[1])/gs[1];
    end if;
    F := Evaluate(G0, mu*x^2);
    Fs := PK ! [ sig(c) : c in Coefficients(F) ];
    Frev := PK ! Reverse(Coefficients(F));          // x^8 F(1/x)
    require Degree(Frev) eq 8 and Degree(Fs) eq 8 : "F must have degree 8 with nonzero constant term";
    cst := Coefficient(Fs, 0)/Coefficient(Frev, 0);
    require Fs eq cst*Frev : "the Galois condition sigma(F)(x) = c x^8 F(1/x) fails";
    // the descent
    PKs<s> := PolynomialRing(K);
    Pt := PKs ! (&+[ Coefficient(F, i) * (1 + w*s)^i * (1 - w*s)^(8-i) : i in [0..8] ]);
    j := Min([ i : i in [0..8] | Coefficient(Pt, i) ne 0 ]);
    Pq := Pt / Coefficient(Pt, j);
    require forall{ c : c in Coefficients(Pq) | c in Rationals() } : "the descended octic is not rational";
    PQ<S> := PolynomialRing(Rationals());
    P := PQ ! [ Rationals() ! c : c in Coefficients(Pq) ];
    dd := LCM([ Denominator(c) : c in Coefficients(P) ]); P := P*dd;
    cc := GCD([ Integers() ! c : c in Coefficients(P) ]); P := P/cc;
    // the quotient: G(t) = G0(mu t) with sigma(t) = 1/t
    Gt := Evaluate(G0, mu*x);
    PEt := PKs ! (&+[ Coefficient(Gt, i) * (1 + w*s)^i * (1 - w*s)^(4-i) : i in [0..4] ]);
    j := Min([ i : i in [0..4] | Coefficient(PEt, i) ne 0 ]);
    PEq := PEt / Coefficient(PEt, j);
    require forall{ c : c in Coefficients(PEq) | c in Rationals() } : "the descended quartic is not rational";
    PE := PQ ! [ Rationals() ! c : c in Coefficients(PEq) ];
    dd := LCM([ Denominator(c) : c in Coefficients(PE) ]); PE := PE*dd;
    cc := GCD([ Integers() ! c : c in Coefficients(PE) ]); PE := PE/cc;
    return P, PE, mu, F;
end intrinsic;

intrinsic NelsonQuadraticTwist(D::RngIntElt, N::RngIntElt, W::SetEnum, um::RngIntElt,
                               A::FldNumElt, B::FldRatElt : PrimeCap := 150)
          -> FldRatElt, FldRatElt, Rec
{The form over Q of the marked curve with K-invariants A, B (NelsonQuadraticDescent) that is
 X = X_0^D(N)/W, found from the trace formula.  Over K the forms are y^2 = kappa G0((mu/e) x^2)
 with the descent datum (x, y) -> (e/x, e^2 y/x^4), where e in Q^x/(Q^x)^2 is the twist by the
 involution u and kappa = r kappa_0 with kappa_0 the Hilbert-90 scalar the datum requires and r
 the quadratic twist.  At a prime split in K the curve over F_p is the K-model reduced at a prime
 above p, and at an inert prime its F_p-points are the points over F_(p^2) fixed by the Frobenius
 composed with the datum.  e is found first, from the split primes, as the class whose counts
 match for some quadratic twist (a linear system over F_2 in its exponents over the primes of
 2DNd and of the norm-equation solution mu); then r from the split primes on y^2 = kappa_0 F_e;
 the factor d in each, invisible at split primes, is settled at the inert primes, where it is the
 sign of the datum.  Returns e, r and a record with: the K-model F = r kappa_0 F_e over K with the
 datum (e, eps), the conic X^2 - d Y^2 = e Z^2 that X/iota is over Q and whether it has a rational
 point (e a norm from K), and when it has, the model y^2 = rQ P(s) over Q (rQ its own twist, in
 the normalisation of P) with its reduced minimal form, the quotient y^2 = r_E P_E(v) over Q, and
 the number of point counts the model matched.  When the conic has no rational point the curve
 has no model y^2 = f(s) over Q and the K-model with its datum is the result.  The counts fix the
 isogeny class of the Jacobian and nothing more.}
    DN := D*N;
    K := Parent(A);
    require Degree(K) eq 2 : "A must lie in a quadratic field";
    w := K.1; d := Integers() ! (w^2);
    sig := hom< K -> K | -w >;
    OK := Integers(K);
    WU := { ALMul(a, um) : a in W } join W;
    require GenusShimuraCurveQuotient(D, N, W) eq 3 and GenusShimuraCurveQuotient(D, N, WU) eq 1 :
        "expected genus 3 over a genus-1 Atkin-Lehner quotient";
    XQ := CreateShimuraQuot(D, N, W);  XQ`g := 3;
    EQ := CreateShimuraQuot(D, N, WU); EQ`g := 1;
    P1, PE, mu, F1 := NelsonQuadraticDescent(A, B);
    PK<x> := PolynomialRing(K);
    G0 := 1 + x + A*x^2 + A^2*B*x^3 + A^3*B^2*sig(A)*x^4;
    // bad primes: those of 2DN, of K, and of the discriminant and coefficients of G0(mu x^2)
    Fmu := Evaluate(G0, mu*x^2);
    cs := Coefficients(Fmu);
    den := LCM([ Denominator(c) : c in cs ]);
    bad := Set(PrimeDivisors(2*DN*Abs(d))) join Set(PrimeDivisors(den))
           join Set(PrimeDivisors(Integers() ! Norm(den^16*Discriminant(Fmu))));
    for c in cs do if c ne 0 then bad join:= Set(PrimeDivisors(Integers() ! Norm(den*c))); end if; end for;
    supp := Sort(SetToSequence(bad));
    nsupp := #supp + 1;   // the sign is the last unknown
    // reduction of K at a split prime
    function Reduce(c, pi, red)
        dn := Denominator(c);
        return red(OK ! (c*dn)) / red(OK ! dn);
    end function;
    // the coefficients and mu must be p-integral for the reduction to mean anything
    function Reducible(p)
        return forall{ c : c in Coefficients(G0) cat [mu] | Denominator(c) mod p ne 0 };
    end function;
    function CountSplit(eps_e, eps_r, p, pi, red, np)
        if not Reducible(p) then return -1; end if;
        k := Codomain(red);
        Pk := PolynomialRing(k);
        eta := eps_e select k ! 1 else k ! np;      // e a square at p or not
        rho := eps_r select k ! 1 else k ! np;
        // y^2 = rho G0((mu/e) x^2) with mu/e reduced: mu/eta
        mred := Reduce(mu, pi, red)/eta;
        g := Pk ! [ Reduce(Coefficient(G0, i), pi, red) : i in [0..4] ];
        f := rho * Evaluate(g, mred * Pk.1^2);
        if Degree(f) lt 7 or Discriminant(f) eq 0 then return -1; end if;
        return #Points(HyperellipticCurve(f));
    end function;
    tX := AssociativeArray(); tE := AssociativeArray();
    good := [ p : p in PrimesUpTo(PrimeCap) | p ne 2 and DN mod p ne 0 ];
    for p in good do
        tX[p] := [ ComputePointsViaTrace(XQ, p, dg) : dg in [1,2] ];
        tE[p] := [ ComputePointsViaTrace(EQ, p, dg) : dg in [1,2] ];
    end for;
    splitp := [ p : p in good | KroneckerSymbol(Discriminant(OK), p) eq 1 ];
    inertp := [ p : p in good | KroneckerSymbol(Discriminant(OK), p) eq -1 ];
    // e is found before r, because the quadratic twist of the K-model is a K-number until the
    // model is normalised by the scalar the descent datum requires.  e is supported on the primes
    // of 2 D N d and on those below the prime ideals dividing (mu): the solution of the norm
    // equation is defined only up to norm-one elements, which shift e by norms.
    mufac := Factorization(mu*OK);
    suppE := Sort(SetToSequence(Set(PrimeDivisors(2*DN*Abs(d))) join { Minimum(f[1]) : f in mufac }));
    nsuppE := #suppE + 1;
    dvec := [ GF(2) | (Abs(d) mod q eq 0 select 1 else 0) : q in suppE ] cat [ GF(2) | (d lt 0 select 1 else 0) ];
    function SymRow(p)
        return [ GF(2) | (LegendreSymbol(q, p) eq 1 select 0 else 1) : q in suppE ] cat [ GF(2) | (p mod 4 eq 1 select 0 else 1) ];
    end function;
    function FromVec(vec)
        sgn := (vec[nsuppE] eq GF(2)!1) select -1 else 1;
        return Rationals() ! (sgn * (&*([Integers() | 1] cat [ suppE[i] : i in [1..#suppE] | vec[i] eq GF(2)!1 ])));
    end function;
    rows := []; rhs := []; nsplit := 0;
    for p in splitp do
        pi := Factorization(p*OK)[1][1];
        k, red := ResidueClassField(pi);
        np := 2; while LegendreSymbol(np, p) eq 1 do np +:= 1; end while;
        adm := { };
        for ee in [true, false], rr in [true, false] do
            c := CountSplit(ee, rr, p, pi, red, np);
            if c eq tX[p][1] then Include(~adm, ee); end if;
        end for;
        if #adm ne 1 then continue; end if;
        nsplit +:= 1;
        Append(~rows, SymRow(p)); Append(~rhs, GF(2) ! (Rep(adm) select 0 else 1));
    end for;
    require nsplit ge 4 : "too few split primes decide the twist by u; raise PrimeCap";
    M := Matrix(GF(2), #rows, nsuppE, &cat rows);
    ok, v := IsConsistent(Transpose(M), Vector(GF(2), rhs));
    require ok : "the split-prime counts are inconsistent with every twist by u: the invariants or the field are wrong";
    Kmat := Kernel(Transpose(M));
    require Dimension(Kmat) le 1 and forall{ k : k in Basis(Kmat) | k eq Vector(GF(2), dvec) } :
        Sprintf("the twist by u is not determined by the split primes beyond the factor d (kernel dimension %o); raise PrimeCap", Dimension(Kmat));
    e0 := FromVec(v);
    vprintf NelsonShimizu, 2: "split primes give e = %o modulo the factor d = %o\n", e0, d;
    // The datum (x, y) -> (e/x, e^2 y/x^4) forces sigma(kappa)/kappa = e^4/c, of norm one, so
    // Hilbert 90 supplies kappa_0 = 1 + c/e^4; every valid scalar is kappa_0 times a rational
    // number (the quadratic twist) times a square.  This needs no e that is a norm.
    function Datum(ee)
        Fe := Evaluate(G0, (mu/ee)*x^2);
        Fs := PK ! [ sig(c) : c in Coefficients(Fe) ];
        num := PK ! (&+[ Coefficient(Fs, i) * ee^i * x^(8-i) : i in [0..8] ]);   // sigma(F)(e/x) x^8
        cst := Coefficient(num, 0)/Coefficient(Fe, 0);
        if num ne cst*Fe then return false, Fe, K ! 0, K ! 0; end if;
        eta := ee^4/cst;
        assert Norm(eta) eq 1;
        kap0 := 1 + 1/eta;
        if kap0 eq 0 then kap0 := w + sig(w)/eta; end if;
        assert sig(kap0)/kap0 eq eta;
        return true, kap0*Fe, K ! ee^2, kap0;
    end function;
    // count over F_p at an inert prime: x in F_{p^2} with x^{p+1} = e, y^2 = F(x), y^p = eps y/x^4
    function CountInert(poly, ee, eps, p)
        if not Reducible(p) or Denominator(ee) mod p eq 0 or Numerator(ee) mod p eq 0 then return -1; end if;
        if exists{ c : c in Coefficients(poly) cat [eps] | Denominator(c) mod p eq 0 } then return -1; end if;
        pid := Factorization(p*OK)[1][1];
        k, red := ResidueClassField(pid);
        function rd(c) return red(OK ! (c*Denominator(c)))/red(OK ! Denominator(c)); end function;
        Pk := PolynomialRing(k);
        f := Pk ! [ rd(Coefficient(poly, i)) : i in [0..8] ];
        if Degree(f) lt 7 or Discriminant(f) eq 0 then return -1; end if;
        epsk := rd(eps); ek := k ! (Integers() ! (GF(p) ! ee));
        cnt := 0;
        for xx in k do
            if xx eq 0 or xx^(p+1) ne ek then continue; end if;
            fx := Evaluate(f, xx);
            if fx eq 0 then cnt +:= 1; continue; end if;
            // y^2 = fx and y^(p-1) = eps/x^4: the second reads fx^((p-1)/2) = eps/x^4, the same
            // for both square roots, so either both count or neither
            if fx^((p-1) div 2) eq epsk/xx^4 then cnt +:= 2; end if;
        end for;
        return cnt;
    end function;
    // count over F_p at a split prime of the normalised model y^2 = rho Fhat
    function CountSplitN(Fhat, rho, p)
        if exists{ c : c in Coefficients(Fhat) | Denominator(c) mod p ne 0 and false } then return -1; end if;
        if exists{ c : c in Coefficients(Fhat) | Denominator(c) mod p eq 0 } then return -1; end if;
        pi := Factorization(p*OK)[1][1]; k, red := ResidueClassField(pi); Pk := PolynomialRing(k);
        f := (k ! rho) * (Pk ! [ Reduce(Coefficient(Fhat, i), pi, red) : i in [0..8] ]);
        if Degree(f) lt 7 or Discriminant(f) eq 0 then return -1; end if;
        return #Points(HyperellipticCurve(f));
    end function;
    cands := [];
    for ee in [e0, Rationals() ! (e0*d)], sg in [1, -1] do
        okd, Fhat, eps, kap0 := Datum(ee);
        vprintf NelsonShimizu, 2: "   candidate e = %o, sign %o: datum exists %o\n", ee, sg, okd;
        if not okd then continue; end if;
        // r is rational relative to kappa_0, whose own rational content is arbitrary, so the
        // primes of N(kappa_0) join the support
        nk := Norm(kap0);
        suppR := Sort(SetToSequence(Set(supp) join Set(PrimeDivisors(Numerator(nk))) join Set(PrimeDivisors(Denominator(nk)))));
        nsuppR := #suppR + 1;
        dvecR := [ GF(2) | (Abs(d) mod q eq 0 select 1 else 0) : q in suppR ] cat [ GF(2) | (d lt 0 select 1 else 0) ];
        rws := []; rh := [];
        for p in splitp do
            n1 := CountSplitN(Fhat, 1, p);
            if n1 lt 0 then continue; end if;
            np := 2; while LegendreSymbol(np, p) eq 1 do np +:= 1; end while;
            n2 := CountSplitN(Fhat, np, p);
            t1 := tX[p][1];
            if n1 eq t1 and n2 eq t1 then continue; end if;
            if n1 eq t1 then b := 0; elif n2 eq t1 then b := 1; else rws := []; break; end if;
            Append(~rws, [ GF(2) | (LegendreSymbol(q, p) eq 1 select 0 else 1) : q in suppR ] cat [ GF(2) | (p mod 4 eq 1 select 0 else 1) ]);
            Append(~rh, GF(2) ! b);
        end for;
        if #rws lt 4 then continue; end if;
        Mr := Matrix(GF(2), #rws, nsuppR, &cat rws);
        okr, vr := IsConsistent(Transpose(Mr), Vector(GF(2), rh));
        if not okr then vprintf NelsonShimizu, 2: "      no rational twist fits the split primes\n"; continue; end if;
        Kr := Kernel(Transpose(Mr));
        if Dimension(Kr) gt 1 or exists{ k : k in Basis(Kr) | k ne Vector(GF(2), dvecR) } then
            vprintf NelsonShimizu, 2: "      the quadratic twist is not determined by the split primes (kernel %o)\n", Dimension(Kr); continue;
        end if;
        sgnr := (vr[nsuppR] eq GF(2)!1) select -1 else 1;
        rr := Rationals() ! (sgnr * (&*([Integers() | 1] cat [ suppR[i] : i in [1..#suppR] | vr[i] eq GF(2)!1 ])));
        agree := true; nin := 0;
        for p in inertp do
            c := CountInert(rr*Fhat, ee, sg*eps, p);
            if c lt 0 then continue; end if;
            nin +:= 1;
            if c ne tX[p][1] then
                vprintf NelsonShimizu, 2: "      inert p = %o: count %o against trace %o\n", p, c, tX[p][1];
                agree := false; break;
            end if;
        end for;
        if agree and nin ge 3 then Append(~cands, <ee, sg*eps, rr, Fhat>); end if;
    end for;
    require #cands ge 1 : "no candidate (e, sign of eps) matches the counts at the inert primes";
    require #cands eq 1 : Sprintf("several candidates match: %o", [ c[1] : c in cands ]);
    e := cands[1][1]; eps := cands[1][2]; r := cands[1][3]; Fhat := cands[1][4];
    vprintf NelsonShimizu, 2: "   chosen e = %o, r = %o (relative to the Hilbert-90 scalar), eps = %o\n", e, r, eps;
    Fe := Fhat;
    con := Conic([ Rationals() | 1, -d, -e ]);
    split := HasRationalPoint(con);
    R := rec< recformat< e, r, rQ, rE, K, F, conic, split, P, PE, model, quotient, counts, mu, eps > |
              e := e, r := r, K := K, F := r*Fe, conic := con, split := split, mu := mu, eps := eps >;
    vprintf NelsonShimizu, 1: "u-twist e = %o, quadratic twist r = %o; conic x^2 - %o y^2 = %o z^2 has a rational point: %o\n", e, r, d, e, split;
    // the quotient by u is the same for every e
    function SolveE()
        suppE := PolySupport(PE, 2*DN*Abs(d));
        rws := []; rh := [];
        for p in good do
            n1 := TwistCountQ(PE, p, 1);
            if n1 lt 0 then continue; end if;
            t1 := tE[p][1];
            if n1 eq t1 and 2*(p+1) - n1 eq t1 then continue; end if;
            if n1 eq t1 then b := 0; elif 2*(p+1) - n1 eq t1 then b := 1; else continue; end if;
            Append(~rws, TwistSymbolRow(suppE, p)); Append(~rh, GF(2) ! b);
        end for;
        return TwistSolve(suppE, rws, rh, "r_E");
    end function;
    rE := SolveE();
    R`rE := rE; R`PE := PE; R`quotient := HyperellipticCurve(rE*PE);
    for p in good do for dg in [1,2] do
        c := TwistCountQ(rE*PE, p, dg);
        require c lt 0 or c eq tE[p][dg] : Sprintf("the quotient's count disagrees at p = %o, degree %o", p, dg);
    end for; end for;
    if split then
        P := NelsonQuadraticDescent(A, B : Twist := e);
        // r exactly, from all odd good primes on the model over Q
        suppR := PolySupport(P, 2*DN*Abs(d));
        rws := []; rh := [];
        for p in good do
            n1 := TwistCountQ(P, p, 1);
            if n1 lt 0 then continue; end if;
            t1 := tX[p][1];
            if n1 eq t1 and 2*(p+1) - n1 eq t1 then continue; end if;
            if n1 eq t1 then b := 0; elif 2*(p+1) - n1 eq t1 then b := 1; else continue; end if;
            Append(~rws, TwistSymbolRow(suppR, p)); Append(~rh, GF(2) ! b);
        end for;
        rQ := TwistSolve(suppR, rws, rh, "r");
        // R`F keeps the K-model with its own twist; rQ is the twist of the model over Q, whose
        // normalisation differs from the K-model's by the constant of the descent
        R`rQ := rQ;
        Cfin := HyperellipticCurve(rQ*P);
        try Cfin := ReducedMinimalWeierstrassModel(Cfin); catch err
            try Cfin := ReducedModel(SimplifiedModel(Cfin)); catch err2 ; end try;
        end try;
        require Genus(Cfin) eq 3 : "genus mismatch";
        fq, hq := HyperellipticPolynomials(Cfin);
        gg := fq + hq^2/4;
        nchk := 0;
        for p in good do for dg in [1,2] do
            c := TwistCountQ(gg, p, dg);
            if c lt 0 then continue; end if;
            require c eq tX[p][dg] : Sprintf("point count disagrees at p = %o, degree %o: %o against %o", p, dg, c, tX[p][dg]);
            nchk +:= 1;
        end for; end for;
        R`P := P; R`model := Cfin; R`counts := nchk;
        vprintf NelsonShimizu, 1: "model over Q: %o (%o point counts agree)\n", Cfin, nchk;
    else
        // no model y^2 = f(s) over Q; count at split primes on the K-model, degrees 1 and 2
        nchk := 0;
        for p in good do
            if KroneckerSymbol(Discriminant(OK), p) ne 1 then continue; end if;
            if not Reducible(p) or Denominator(r) mod p eq 0 or Denominator(e) mod p eq 0 or Numerator(e) mod p eq 0 then continue; end if;
            pi := Factorization(p*OK)[1][1]; k, red := ResidueClassField(pi);
            Pk := PolynomialRing(k);
            f := Pk ! [ Reduce(Coefficient(r*Fe, i), pi, red) : i in [0..8] ];
            if Degree(f) lt 7 or Discriminant(f) eq 0 then continue; end if;
            C := HyperellipticCurve(f);
            for dg in [1,2] do
                c := #Points(BaseChange(C, GF(p^dg)));
                require c eq tX[p][dg] : Sprintf("K-model count disagrees at p = %o, degree %o", p, dg);
                nchk +:= 1;
            end for;
        end for;
        for p in good do
            if KroneckerSymbol(Discriminant(OK), p) ne -1 then continue; end if;
            c := CountInert(r*Fe, e, eps, p);
            if c lt 0 then continue; end if;
            require c eq tX[p][1] : Sprintf("K-model count disagrees at the inert prime %o", p);
            nchk +:= 1;
        end for;
        R`counts := nchk;
        vprintf NelsonShimizu, 1: "no model y^2 = f(s) over Q: the conic is pointless; K-model y^2 = %o with sigma(x) = %o/x (%o counts agree)\n", r*Fe, e, nchk;
    end if;
    return e, r, R;
end intrinsic;

intrinsic NelsonRecogniseInQuadraticField(A::FldComElt, C::FldComElt, d::RngIntElt, eps::FldReElt) -> BoolElt, FldNumElt, RngIntElt
{Recognise A as an element (p + q sqrt d)/r of K = Q(sqrt d), with C its conjugate (for d < 0 the
 complex conjugate of A up to noise; for d > 0 a second real number), from the two rationals
 (A + C)/2 and (A - C)/(2 sqrt d) at absolute accuracy eps, by lattice reduction on the pair.
 A solution is accepted when both residuals are below eps r and r is below the bound
 (1/(2 eps^2))^(1/3), the denominator up to which a spurious simultaneous approximation of two
 numbers to that accuracy is not expected; the bound is returned as the third value.  As with
 any such recognition, the result is certain only relative to that bound, and an element of K of
 larger height is missed, not misread.}
    RR := Parent(eps);
    K := QuadraticField(d); w := K.1;
    sd := Sqrt(RR ! Abs(d));
    if d lt 0 then
        re := RR ! (Real(A) + Real(C))/2;
        im := RR ! (Imaginary(A) - Imaginary(C))/(2*sd);
    else
        re := RR ! (Real(A) + Real(C))/2;
        im := RR ! (Real(A) - Real(C))/(2*sd);
    end if;
    bound := Floor((1/(2*eps^2))^(1/3));
    S := Round(RR ! 1/eps);
    L := LatticeWithBasis(Matrix(Integers(), 3, 5, [
        1, 0, 0, Round(S*re), Round(S*im),
        0, 1, 0, -S, 0,
        0, 0, 1, 0, -S ]));
    Lr := LLL(L);
    for i in [1..3] do
        v := Lr.i;
        r := Integers() ! v[1]; p := Integers() ! v[2]; q := Integers() ! v[3];
        if r eq 0 then continue; end if;
        if r lt 0 then r := -r; p := -p; q := -q; end if;
        if r gt bound then continue; end if;
        if Abs(r*re - p) lt eps*r and Abs(r*im - q) lt eps*r then
            return true, (p + q*w)/r, bound;
        end if;
    end for;
    return false, K ! 0, bound;
end intrinsic;

intrinsic NelsonBranchField(D::RngIntElt, N::RngIntElt, W::SetEnum, um::RngIntElt) -> RngIntElt, SeqEnum
{The discriminants of the quadratic subfields of the fields of definition of the fixed points of
 u = w_um on X = X_0^D(N)/W (the CM points of the coset u W), and the list <discriminant, points
 on X/W, degree of the field, quadratic subfield discriminants>.  When the fixed points of u on
 the canonical conic are conjugate over K, every fixed point of u on X is defined over a field
 containing K, so K is a common quadratic subfield; the first return value is the discriminant
 of the unique common quadratic subfield, or 0 when the fields have degree at most 2 with
 different or no quadratic subfield (the rational case).}
    DN := D*N;
    WU := { ALMul(a, um) : a in W } join W;
    coset := WU diff W;
    XQ := CreateShimuraQuot(D, N, W); XQ`g := 3;
    out := []; common := {};
    first := true;
    for m in Sort(SetToSequence(coset)) do
        fp := NumFixedPointsByCMOrder(D, N, m);
        for dsc in Sort(SetToSequence(Keys(fp))) do
            if fp[dsc] eq 0 then continue; end if;
            F := FieldOfDefinitionOfCMPoint(XQ, dsc);
            deg := Type(F) eq FldRat select 1 else Degree(F);
            subs := deg eq 1 select [] else (deg eq 2 select [ Discriminant(Integers(F)) ]
                    else [ Discriminant(Integers(s[1])) : s in Subfields(F) | Degree(s[1]) eq 2 ]);
            Append(~out, <dsc, fp[dsc]/#W, deg, subs>);
            if first then common := Set(subs); first := false; else common meet:= Set(subs); end if;
        end for;
    end for;
    degs := [ o[3] : o in out ];
    if #common eq 1 and Max(degs) eq 4 then return Rep(common), out; end if;
    if #common eq 1 and Max(degs) eq 2 and #out eq 1 then return Rep(common), out; end if;
    return 0, out;
end intrinsic;

/////////////////////////////////////////////////////////////////////////////////////////////////
// The checks beyond point counts

intrinsic NelsonDescentChecks(D::RngIntElt, N::RngIntElt, W::SetEnum, um::RngIntElt, R::Rec :
                              ZetaCap := 60, PrimeCap := 100) -> BoolElt, Rec
{The checks (a)-(d) of NelsonModelChecks for a model found by NelsonQuadraticTwist (R its record).
 (a) The four fixed points of u on X lie over x = 0 and x = infinity of the K-coordinate; on the
 K-model y^2 = F(x) those over x = 0 are defined over K(sqrt F(0)), a field of degree 2 or 4
 over Q (the ones over infinity over its conjugate), which must be isomorphic to the field of
 definition of the CM points of the coset u W, and their total over #W must be 4.  (b), (c) as
 before, with the Jacobian of the descended quotient y^2 = r_E P_E.  (d) The numerator of the zeta
 function over F_p from the K-model reduced at a prime above a SPLIT p <= ZetaCap against the trace
 formula's counts over F_p, F_p^2, F_p^3 (at inert primes only the count over F_p is compared,
 through the descent datum); and the quotient's likewise at every good p.  Numerical plus point
 counts, as the model itself.}
    DN := D*N;
    K := R`K; w := K.1; d := Integers() ! (w^2); OK := Integers(K); sig := hom< K -> K | -w >;
    F := R`F; e := R`e;
    WU := { ALMul(a, um) : a in W } join W;
    coset := WU diff W;
    out := rec< recformat< branchmodel, branchcm, splitting, jacobian, newformap, jacobianap, zetachecked, zetafailed, failures > | failures := [] >;
    fails := [];
    // (a)
    PK<x> := PolynomialRing(K);
    f0 := Coefficient(F, 0);
    issq, rt := IsSquare(f0);
    Fmodel := issq select K else AbsoluteField(ext< K | Polynomial(K, [-f0, 0, 1]) >);
    out`branchmodel := Fmodel;
    tot := 0; cmfields := [* *];
    XQ := CreateShimuraQuot(D, N, W); XQ`g := 3;
    for m in Sort(SetToSequence(coset)) do
        tot +:= NumFixedPoints(D, N, m);
        fp := NumFixedPointsByCMOrder(D, N, m);
        for dsc in Sort(SetToSequence(Keys(fp))) do
            if fp[dsc] gt 0 then Append(~cmfields, <dsc, fp[dsc]/#W, FieldOfDefinitionOfCMPoint(XQ, dsc)>); end if;
        end for;
    end for;
    out`branchcm := cmfields;
    if tot ne 4*#W then Append(~fails, Sprintf("(a) %o fixed points of the coset over #W = %o, not 4", tot, #W)); end if;
    for c in cmfields do
        Fc := c[3];
        same := (Type(Fc) eq FldRat and Degree(Fmodel) eq 1) or
                (Type(Fc) ne FldRat and Degree(Fc) eq Degree(Fmodel) and IsIsomorphic(Fc, Fmodel));
        if not same then Append(~fails, Sprintf("(a) the CM points of discriminant %o are defined over %o, the model's fixed points over %o", c[1], Fc, Fmodel)); end if;
    end for;
    // (b), (c)
    ps := [ p : p in PrimesUpTo(PrimeCap) | DN mod p ne 0 ];
    plusdim := 0; minusdim := 0; plusf := [];
    for orb in Newforms(CuspForms(Gamma0(DN), 2)) do
        f := orb[1];
        if Level(Parent(f)) ne DN then continue; end if;
        dg := Degree(BaseRing(Parent(f)));
        q := qExpansion(f, Max(PrimeDivisors(DN)) + 2);
        ok := true; eu := 0;
        for m in Sort(SetToSequence(WU)) do
            if m eq 1 then continue; end if;
            ev := 1;
            for p in PrimeDivisors(m) do ev *:= (D mod p eq 0) select Coefficient(q, p) else -Coefficient(q, p); end for;
            if m in W and ev ne 1 then ok := false; end if;
            if m eq um then eu := ev; end if;
        end for;
        if not ok then continue; end if;
        if eu eq 1 then plusdim +:= dg; Append(~plusf, f); else minusdim +:= dg; end if;
    end for;
    out`splitting := [plusdim, minusdim];
    if not (plusdim eq 1 and minusdim eq 2) then Append(~fails, Sprintf("(c) splitting %o + %o", plusdim, minusdim)); end if;
    E := MinimalModel(Jacobian(GenusOneModel(R`rE*R`PE)));
    out`jacobian := E;
    if Conductor(E) ne DN then Append(~fails, Sprintf("(b) conductor %o", Conductor(E))); end if;
    if #plusf eq 1 then
        q := qExpansion(plusf[1], PrimeCap + 1);
        af := [ Integers() ! Coefficient(q, p) : p in ps ];
        aE := [ TraceOfFrobenius(E, p) : p in ps ];
        out`newformap := af; out`jacobianap := aE;
        if af ne aE then Append(~fails, "(b) the Jacobian's a_p differ from the newform's"); end if;
    end if;
    // (d)
    EQ := CreateShimuraQuot(D, N, WU); EQ`g := 1;
    RT<T> := PolynomialRing(Rationals());
    function TraceL(Xq, g, p)
        cnt := [ ComputePointsViaTrace(Xq, p, r) : r in [1..g] ];
        s := [ p^r + 1 - cnt[r] : r in [1..g] ];
        ee := [ Rationals() | ];
        for r in [1..g] do
            acc := s[r];
            for j in [1..r-1] do acc +:= (-1)^j * ee[j] * s[r-j]; end for;
            Append(~ee, (-1)^(r-1)*acc/r);
        end for;
        a := [ Rationals() | 1 ] cat [ (-1)^r*ee[r] : r in [1..g] ];
        for r in [g-1..0 by -1] do Append(~a, p^(g-r)*a[r+1]); end for;
        return &+[ a[i+1]*T^i : i in [0..2*g] ];
    end function;
    nchk := 0; nfz := 0;
    for p in [ p : p in PrimesUpTo(ZetaCap) | (2*DN*Abs(d)) mod p ne 0 ] do
        if KroneckerSymbol(Discriminant(OK), p) ne 1 then continue; end if;
        if exists{ c : c in Coefficients(F) | Denominator(c) mod p eq 0 } then continue; end if;
        pid := Factorization(p*OK)[1][1]; k, red := ResidueClassField(pid); Pk := PolynomialRing(k);
        fp := Pk ! [ red(OK ! (c*Denominator(c)))/red(OK ! Denominator(c)) : c in Coefficients(F) ];
        if Degree(fp) lt 7 or Discriminant(fp) eq 0 then continue; end if;
        nchk +:= 1;
        if RT ! LPolynomial(HyperellipticCurve(fp)) ne TraceL(XQ, 3, p) then nfz +:= 1; end if;
    end for;
    for p in [ p : p in PrimesUpTo(ZetaCap) | (2*DN) mod p ne 0 ] do
        if TwistCountQ(R`rE*R`PE, p, 1) lt 0 then continue; end if;
        nchk +:= 1;
        if 1 - TraceOfFrobenius(E, p)*T + p*T^2 ne TraceL(EQ, 1, p) then nfz +:= 1; end if;
    end for;
    out`zetachecked := nchk; out`zetafailed := nfz;
    if nfz gt 0 then Append(~fails, Sprintf("(d) %o of %o zeta numerators differ", nfz, nchk)); end if;
    out`failures := fails;
    vprintf NelsonShimizu, 1: "checks: fixed points over %o; CM %o; splitting %o; Jacobian %o; zeta %o/%o; %o\n",
        Fmodel, [ <c[1], c[2], c[3]> : c in cmfields ], out`splitting, E, nchk - nfz, nchk, #fails eq 0 select "all pass" else fails;
    return #fails eq 0, out;
end intrinsic;

intrinsic NelsonModelChecks(D::RngIntElt, N::RngIntElt, W::SetEnum, um::RngIntElt,
                            ABC::SeqEnum, twist::SeqEnum : ZetaCap := 60, PrimeCap := 100) -> BoolElt, Rec
{Four checks of a model y^2 = kappa g0(mu x^2) of X = X_0^D(N)/W produced by
 NelsonHyperellipticModel, given the invariants ABC = [A, B, C] and twist = [kappa, mu], with
 u = w_um acting as x -> -x and X/<W,u> : Y^2 = kappa g0(t).
 (a) The four fixed points of u on X are the images of the fixed points of the elements of the
 coset uW (Atkin-Lehner group law), CM points whose discriminants and counts come from
 NumFixedPointsByCMOrder; their total over #W must be 4, and the fields over which they are
 defined on X (FieldOfDefinitionOfCMPoint) must match the model, where the pair over x = 0 is
 defined over Q(sqrt kappa) and the pair over x = infinity over Q(sqrt of the leading
 coefficient): each CM discriminant contributes one entry per pair, so a discriminant carrying
 four points of the quotient defined over one quadratic field contributes that field twice.
 (b) The Jacobian of X/<W,u> (Magma's Jacobian of the quartic) has conductor DN and the same a_p,
 p < PrimeCap, as the W-fixed newform of level DN with quaternionic u-eigenvalue +1.
 (c) The W-fixed newforms of level DN split 1 + 2 under u.
 (d) The numerator of the zeta function of the model over F_p, for good p up to ZetaCap, equals
 the one assembled from the trace formula's counts over F_p, F_p^2, F_p^3; likewise for the
 quotient.  Returns whether all pass and a record of the findings.  These determine the isogeny
 class of the Jacobian and the Galois structure of the branch locus, not the curve.}
    require #ABC eq 3 and #twist eq 2 : "three invariants and two twist constants expected";
    A := Rationals() ! ABC[1]; B := Rationals() ! ABC[2]; C := Rationals() ! ABC[3];
    kappa := Rationals() ! twist[1]; mu := Rationals() ! twist[2];
    PQ<x> := PolynomialRing(Rationals());
    DN := D*N;
    gnorm := 1 + x + A*x^2 + A^2*B*x^3 + A^3*B^2*C*x^4;
    foct := kappa*Evaluate(gnorm, mu*x^2);
    fE := kappa*gnorm;
    WU := { ALMul(a, um) : a in W } join W;
    coset := WU diff W;
    R := rec< recformat< branchmodel, branchcm, branchdiscs, splitting, jacobian, newformap, jacobianap,
                         zetachecked, zetafailed, failures > | failures := [] >;
    fails := [];

    // (a)
    sfk := Squarefree(Numerator(kappa)*Denominator(kappa));
    lead := LeadingCoefficient(foct);
    sfl := Squarefree(Numerator(lead)*Denominator(lead));
    predicted := {* sfk, sfl *};
    R`branchmodel := predicted;
    tot := 0; discs := [];
    for m in Sort(SetToSequence(coset)) do
        tot +:= NumFixedPoints(D, N, m);
        fp := NumFixedPointsByCMOrder(D, N, m);
        for d in Sort(SetToSequence(Keys(fp))) do
            if fp[d] gt 0 then Append(~discs, <m, d, fp[d]>); end if;
        end for;
    end for;
    R`branchdiscs := discs;
    if tot ne 4*#W then Append(~fails, Sprintf("(a) %o fixed points of the coset over #W = %o, not 4", tot, #W)); end if;
    XQ := CreateShimuraQuot(D, N, W); XQ`g := 3;
    found := {* *}; badfield := false;
    for dc in discs do
        d := dc[2];
        npts := dc[3] / #W;
        F := FieldOfDefinitionOfCMPoint(XQ, d);
        deg := Type(F) eq FldRat select 1 else Degree(F);
        if deg gt 2 or not IsIntegral(npts) or (Integers() ! npts) mod 2 ne 0 then badfield := true; continue; end if;
        s := deg eq 1 select 1 else Squarefree(Discriminant(F));
        for i in [1..(Integers() ! npts) div 2] do Include(~found, s); end for;
    end for;
    R`branchcm := found;
    if badfield or found ne predicted then
        Append(~fails, Sprintf("(a) branch fields: CM theory %o, model %o", found, predicted));
    end if;

    // (b), (c)
    ps := [ p : p in PrimesUpTo(PrimeCap) | DN mod p ne 0 ];
    plusdim := 0; minusdim := 0; plusf := [];
    for orb in Newforms(CuspForms(Gamma0(DN), 2)) do
        f := orb[1];
        if Level(Parent(f)) ne DN then continue; end if;
        dg := Degree(BaseRing(Parent(f)));
        q := qExpansion(f, Max(PrimeDivisors(DN)) + 2);
        ok := true; eu := 0;
        for m in Sort(SetToSequence(WU)) do
            if m eq 1 then continue; end if;
            e := 1;
            for p in PrimeDivisors(m) do e *:= (D mod p eq 0) select Coefficient(q, p) else -Coefficient(q, p); end for;
            if m in W and e ne 1 then ok := false; end if;
            if m eq um then eu := e; end if;
        end for;
        if not ok then continue; end if;
        if eu eq 1 then plusdim +:= dg; Append(~plusf, f); else minusdim +:= dg; end if;
    end for;
    R`splitting := [plusdim, minusdim];
    if not (plusdim eq 1 and minusdim eq 2) then Append(~fails, Sprintf("(c) splitting %o + %o", plusdim, minusdim)); end if;
    E := MinimalModel(Jacobian(GenusOneModel(fE)));
    R`jacobian := E;
    if Conductor(E) ne DN then Append(~fails, Sprintf("(b) conductor %o", Conductor(E))); end if;
    if #plusf eq 1 then
        q := qExpansion(plusf[1], PrimeCap + 1);
        af := [ Integers() ! Coefficient(q, p) : p in ps ];
        aE := [ TraceOfFrobenius(E, p) : p in ps ];
        R`newformap := af; R`jacobianap := aE;
        if af ne aE then Append(~fails, "(b) the Jacobian's a_p differ from the newform's"); end if;
    end if;

    // (d)
    EQ := CreateShimuraQuot(D, N, WU); EQ`g := 1;
    RT<T> := PolynomialRing(Rationals());
    nchk := 0; nfz := 0;
    for pair in [ <foct, XQ, 3>, <fE, EQ, 1> ] do
        poly := pair[1]; Xq := pair[2]; g := pair[3];
        den := LCM([ Denominator(c) : c in Coefficients(poly) ]);
        Cq := HyperellipticCurve(den^2*poly);
        for p in [ p : p in PrimesUpTo(ZetaCap) | (2*DN*den) mod p ne 0 ] do
            Cp := ChangeRing(Cq, GF(p));
            if not IsNonsingular(Cp) then continue; end if;
            // the genus-1 quartic need not have a rational point; its Jacobian has the same L-polynomial
            Lmodel := g eq 1 select 1 - TraceOfFrobenius(E, p)*T + p*T^2 else RT ! LPolynomial(Cp);
            cnt := [ ComputePointsViaTrace(Xq, p, r) : r in [1..g] ];
            s := [ p^r + 1 - cnt[r] : r in [1..g] ];
            e := [ Rationals() | ];
            for r in [1..g] do
                acc := s[r];
                for j in [1..r-1] do acc +:= (-1)^j * e[j] * s[r-j]; end for;
                Append(~e, (-1)^(r-1)*acc/r);
            end for;
            a := [ Rationals() | 1 ] cat [ (-1)^r*e[r] : r in [1..g] ];
            for r in [g-1..0 by -1] do Append(~a, p^(g-r)*a[r+1]); end for;   // a_{2g-r} = p^{g-r} a_r
            Ltrace := &+[ a[i+1]*T^i : i in [0..2*g] ];
            nchk +:= 1;
            if Ltrace ne Lmodel then nfz +:= 1; end if;
        end for;
    end for;
    R`zetachecked := nchk; R`zetafailed := nfz;
    if nfz gt 0 then Append(~fails, Sprintf("(d) %o of %o zeta numerators differ", nfz, nchk)); end if;
    R`failures := fails;
    vprintf NelsonShimizu, 1: "checks: branch fields %o against %o; splitting %o; Jacobian %o; zeta %o/%o; %o\n",
        found, predicted, R`splitting, E, nchk - nfz, nchk, #fails eq 0 select "all pass" else fails;
    return #fails eq 0, R;
end intrinsic;
