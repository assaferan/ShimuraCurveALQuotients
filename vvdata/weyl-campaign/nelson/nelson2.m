// Nelson, "Evaluating modular forms on Shimura curves", Math. Comp. 84 (2015) 2471-2503
// (arXiv:1210.1243), Theorem 3.1.  Working version: same formula as nelson.m, reorganised so that
// the per-lattice-point work is four complex multiplications and one Chebyshev evaluation.
//
//   (y1 y2)^{k/2} conj(f_B(z1)) f_B(z2)
//      = sum_{d | DN} c_f(d) sum_{alpha in R^(d), nrd(alpha) > 0} a_f(d nrd alpha) e^{ik arg X} W_k(4 pi |X|)
//
// with X = X(sigma_{z1}^{-1} iota(alpha) sigma_{z2}), c_f(d) = mu((d,N))/(d a_f(d)), and a_f the
// CLASSICAL coefficients of the Jacquet-Langlands correspondent of f_B on Gamma_0(DN).
//
// Two things make the inner loop cheap.  X is LINEAR in alpha, so the four numbers X(b_i) for a
// basis determine every X by a dot product and no 2x2 matrices are formed; and the reduced norm
// is the integral quadratic form v -> v G v^t for G[i][j] = trd(b_i conj(b_j))/2, so no quaternion
// multiplication happens inside the loop either.

load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/nelson/wk.m";

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

function TraceGram2(bas, A)
    n := #bas;
    return Matrix(Rationals(), n, n, [ Trace(A ! (bas[i]*bas[j])) : i, j in [1..n] ]);
end function;

function NormGram(bas, A)
    n := #bas;
    return Matrix(Rationals(), n, n,
        [ Trace(A ! (bas[i] * Conjugate(A ! bas[j]))) / 2 : i, j in [1..n] ]);
end function;

// R^(d): locally R_p^* at p | d, R_p elsewhere.  R^*/R has exponent dividing DN, so this is
// R + (DN/d) R^*, and R^* is spanned by the rows of the inverse of the reduced-trace Gram.
function RdBasis2(T, DN, d)
    n := Nrows(T);
    gens := VerticalJoin(IdentityMatrix(Rationals(), n), (DN div d) * T^(-1));
    den := LCM([ Denominator(x) : x in Eltseq(gens) ]);
    return ChangeRing(BasisMatrix(Lattice(ChangeRing(den*gens, Integers()))), Rationals()) / den;
end function;

function SigmaZ2(z, CC)
    x := Real(z); y := Imaginary(z); r := Sqrt(y);
    return Matrix(CC, 2, 2, [ r, x/r, 0, 1/r ]);
end function;
function SigmaZinv2(z, CC)
    x := Real(z); y := Imaginary(z); r := Sqrt(y);
    return Matrix(CC, 2, 2, [ 1/r, -x/r, 0, r ]);
end function;

NelRec := recformat< D, N, k, A, O, mrep, DN, ds, latb, latQ, prec, Wtab, cutmax >;

// Everything that does not depend on the pair of points.  af is the function n -> a_f(n), needed
// here only at the divisors of DN, for the coefficients c_f(d).
function NelsonSetup(D, N, k, prec, cutmax : OO := 0, MREP := 0)
    if Type(OO) eq RngIntElt then
        A := QuaternionAlgebra(D); O := QuaternionOrder(A, N);
    else
        O := OO; A := Algebra(O);
    end if;
    assert Discriminant(O) eq D*N;
    RR := RealField(prec);
    if Type(MREP) eq RngIntElt then
        mrep := FuchsianMatrixRepresentation(A : Precision := prec);
    else
        mrep := MREP;
    end if;
    bas := [ A ! b : b in ZBasis(O) ];
    T := TraceGram2(bas, A);
    DN := D*N; ds := Divisors(DN);
    latb := AssociativeArray(); latQ := AssociativeArray();
    for d in ds do
        B := RdBasis2(T, DN, d);
        bd := [ &+[ B[i][j]*bas[j] : j in [1..4] ] : i in [1..4] ];
        latb[d] := bd;
        latQ[d] := NormGram(bd, A);          // reduced norm as a quadratic form on coordinates
    end for;
    Wk := func< x | WkDirect(k, x, RR) >;
    Wtab := WkPanels(Wk, RR ! 0.4, RR ! (cutmax + 1), RR : ratio := 1.5, deg := prec + 14);
    return rec< NelRec | D := D, N := N, k := k, A := A, O := O, mrep := mrep, DN := DN,
                         ds := ds, latb := latb, latQ := latQ, prec := prec,
                         Wtab := Wtab, cutmax := cutmax >;
end function;

// Theorem 3.1 at (z1, z2), truncated at 4 pi |X| <= cut.  Returns the value, the number of
// lattice points used, the largest n at which a_f was needed, and the largest single term --
// the ratio of that to the value is how many digits of working precision the cancellation eats.
function NelsonEval(dat, z1, z2, af, cut, CC)
    k := dat`k; RR := RealField(Precision(CC)); pi := Pi(RR); i_ := CC.1;
    assert IsEven(k) and cut le dat`cutmax;
    s1i := SigmaZinv2(z1, CC); s2 := SigmaZ2(z2, CC);
    C := cut/(4*pi); Pb := 2*C^2;
    total := CC ! 0; nterms := 0; nmaxused := 0; biggest := RR ! 0;
    SCL := 10^(Min(16, Precision(CC)));
    for d in dat`ds do
        bs := dat`latb[d]; Qn := dat`latQ[d];
        Ms := [ s1i * ChangeRing(dat`mrep(b), CC) * s2 : b in bs ];
        // X is linear in alpha
        Xs := [ ((M[1][1] + M[2][2]) + i_*(M[1][2] - M[2][1]))/2 : M in Ms ];
        // Gram of P(sigma1^{-1} iota(x) sigma2) = (1/2) * Frobenius^2
        G := Matrix(RR, 4, 4, [ RR ! 0 : j in [1..16] ]);
        for a in [1..4] do for b in [a..4] do
            v := (&+[ Real(Ms[a][r][c])*Real(Ms[b][r][c]) : r, c in [1..2] ])/2;
            G[a][b] := v; G[b][a] := v;
        end for; end for;
        GI := Matrix(Integers(), 4, 4, [ Round(SCL*G[a][b]) : a, b in [1..4] ]);
        for a in [1..4] do for b in [a+1..4] do GI[b][a] := GI[a][b]; end for; end for;
        for a in [1..4] do GI[a][a] +:= 8; end for;
        SV := ShortVectors(LatticeWithGram(GI), Ceiling(SCL*Pb*1.000001));
        cf := (CC ! MoebiusMu(GCD(d, dat`N))) / (d * (CC ! af(d)));
        sub := CC ! 0;
        for sv in SV do
            v := Eltseq(sv[1]);
            t := &+[ v[a]*v[b]*Qn[a][b] : a, b in [1..4] ];
            if t le 0 then continue; end if;
            n := d*t;
            assert IsIntegral(n);
            n := Integers() ! n;
            X := v[1]*Xs[1] + v[2]*Xs[2] + v[3]*Xs[3] + v[4]*Xs[4];
            aX := Abs(X);
            x := 4*pi*aX;
            if x gt cut then continue; end if;
            u := X/aX;
            ph := u^k;
            tm := (CC ! af(n)) * ph * WkFast(dat`Wtab, x);
            sub +:= tm;
            bb := 2*Abs(cf)*Abs(tm);
            if bb gt biggest then biggest := bb; end if;
            nterms +:= 1;
            if n gt nmaxused then nmaxused := n; end if;
        end for;
        total +:= 2*cf*sub;                 // alpha and -alpha, k even
    end for;
    return total, nterms, nmaxused, biggest;
end function;
