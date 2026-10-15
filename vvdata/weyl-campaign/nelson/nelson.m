// Nelson, "Evaluating modular forms on Shimura curves", Math. Comp. 84 (2015) 2471-2503
// (arXiv:1210.1243).  Theorem 3.1:
//
//   (y_1 y_2)^{k/2} conj(f_B(z_1)) f_B(z_2)
//        = sum_{n >= 1} a_f(n) sum_{d | d_B N} c_f(d) V_k( sigma_{z1}^{-1} R^{(d)} sigma_{z2} ; n/d )
//
// with  c_f(d) = mu((d,N)) / (d a_f(d)),
//       V_k(L,t) = sum_{alpha in L, det alpha = t} e^{i k arg X(alpha)} W_k(4 pi |X(alpha)|),
//       X(alpha) = (a + d + i(b - c))/2  for alpha = [a b; c d],
//       P(alpha) = (a^2+b^2+c^2+d^2)/2,
//       W_k(x)  = 2^{1-k} sum_{n>=1} n^k ( n x K_{k-1}(n x) - K_k(n x) ),
//       a_f(n)  the CLASSICAL coefficients of the Jacquet-Langlands lift, f = sum a_f(n) e(nz),
//               a_f(1) = 1.
//
// The identity  2|X(alpha)|^2 = P(alpha) + det(alpha)  is what turns the exponential decay of W_k
// into a finite enumeration: every term with |X| <= C has P <= 2C^2, a positive definite condition.

// ---------------------------------------------------------------- W_k

function WkFun(k, x, RR)
    // 2^{1-k} sum_n n^k (n x K_{k-1}(nx) - K_k(nx)); terms fall off like e^{-n x} so the sum is
    // short for the arguments that occur (x >= 4 pi |X|_min, never small in practice).
    eps := RR ! 10^(-(Precision(RR)+5));
    s := RR ! 0;
    n := 1;
    while true do
        u := n * x;
        t := (RR ! n)^k * (u * KBessel(k-1, u) - KBessel(k, u));
        s +:= t;
        if Abs(t) lt eps * (1 + Abs(s)) and n ge 2 then break; end if;
        n +:= 1;
        if n gt 400 then break; end if;
    end while;
    return s / 2^(k-1);
end function;

// ---------------------------------------------------------------- the lattices R^{(d)}

// Z-basis of O as a matrix, and the trace-form Gram  T[i][j] = trd(b_i b_j).
function TraceGram(bas, A)
    n := #bas;
    return Matrix(Rationals(), n, n, [ Trace(A ! (bas[i]*bas[j])) : i, j in [1..n] ]);
end function;

// R^{(d)}: locally R_p^* at p | d and R_p elsewhere.  Since R^*/R has exponent dividing DN,
// R^{(d)} = R + (DN/d) R^*, and R^* has basis the rows of T^{-1} in the coordinates of bas.
// Returned as a 4 x 4 rational matrix whose rows are a Z-basis, in the coordinates of bas.
function RdBasis(T, DN, d)
    n := Nrows(T);
    Ti := T^(-1);
    gens := VerticalJoin(IdentityMatrix(Rationals(), n), (DN div d) * Ti);
    den := LCM([ Denominator(x) : x in Eltseq(gens) ]);
    Lz := Lattice(ChangeRing(den * gens, Integers()));
    B := ChangeRing(BasisMatrix(Lz), Rationals()) / den;
    return B;
end function;

// ---------------------------------------------------------------- the evaluation

// sigma_z = [y^{1/2}, x y^{-1/2}; 0, y^{-1/2}],  sigma_z (i) = z.
function SigmaZ(z, CC)
    x := Real(z); y := Imaginary(z);
    r := Sqrt(y);
    return Matrix(CC, 2, 2, [ r, x/r, 0, 1/r ]);
end function;

function SigmaZinv(z, CC)
    x := Real(z); y := Imaginary(z);
    r := Sqrt(y);
    return Matrix(CC, 2, 2, [ 1/r, -x/r, 0, r ]);
end function;

// Everything needed at a base, computed once.
//   D, N  : the Shimura curve X_0^D(N);  A, O, mrep the algebra, Eichler order and splitting
//   af    : a_f(n) for 1 <= n <= nmax  (classical coefficients, a_f(1) = 1)
function NelsonData(D, N, k, prec : OO := 0, MREP := 0)
    if Type(OO) eq RngIntElt then
        A := QuaternionAlgebra(D);
        O := QuaternionOrder(A, N);
    else
        O := OO; A := Algebra(O);
    end if;
    assert Discriminant(O) eq D*N;
    if Type(MREP) eq RngIntElt then
        mrep := FuchsianMatrixRepresentation(A : Precision := prec);
    else
        mrep := MREP;
    end if;
    bas := [ A ! b : b in ZBasis(O) ];
    T := TraceGram(bas, A);
    DN := D*N;
    ds := Divisors(DN);
    lat := AssociativeArray();
    for d in ds do
        B := RdBasis(T, DN, d);
        lat[d] := [ &+[ B[i][j]*bas[j] : j in [1..4] ] : i in [1..4] ];
    end for;
    return rec< recformat< D, N, k, A, O, mrep, bas, T, DN, ds, lat, prec >
              | D := D, N := N, k := k, A := A, O := O, mrep := mrep, bas := bas,
                T := T, DN := DN, ds := ds, lat := lat, prec := prec >;
end function;

// a_f(n) for n <= nmax, from the classical newform of level DN, weight k, new at every p | D.
// Returns the list and the Atkin-Lehner-type values a_f(d) for d | DN.
function NewformCoeffs(D, N, k, nmax)
    M := CuspForms(Gamma0(D*N), k);
    NF := Newforms(M);
    good := [];
    for orb in NF do
        f := orb[1];
        L := Level(Parent(f)) eq D*N select D*N else Level(Parent(f));
        if L ne D*N then continue; end if;
        // must be new at every p | D (automatic when the newform has level exactly D*N)
        Append(~good, f);
    end for;
    return good;
end function;

// The right-hand side of Theorem 3.1 at (z1,z2).  cut is the cutoff on 4 pi |X|.
// af : function n -> a_f(n) (as a complex number); afd : function d -> a_f(d).
function NelsonPair(dat, z1, z2, af, cut, CC : Verbose := false)
    k := dat`k; RR := RealField(Precision(CC)); pi := Pi(RR);
    assert IsEven(k);
    s1i := SigmaZinv(z1, CC); s2 := SigmaZ(z2, CC);
    C := cut / (4*pi);                  // |X| <= C
    Pbound := 2*C^2;                    // P(alpha) <= 2|X|^2 - det <= 2C^2
    total := CC ! 0;  nterms := 0;
    for d in dat`ds do
        bs := dat`lat[d];
        Ms := [ s1i * ChangeRing(dat`mrep(b), CC) * s2 : b in bs ];
        // Gram of  Q(x) = (1/2) * Frobenius^2 ( sigma1^{-1} iota(x) sigma2 )
        G := Matrix(RR, 4, 4, [ RR ! 0 : i in [1..16] ]);
        for i in [1..4] do for j in [i..4] do
            v := (&+[ Real(Ms[i][r][c]) * Real(Ms[j][r][c]) : r, c in [1..2] ]) / 2;
            G[i][j] := v; G[j][i] := v;
        end for; end for;
        SCL := 10^12;
        GI := Matrix(Integers(), 4, 4, [ Round(SCL * G[i][j]) : i, j in [1..4] ]);
        // symmetrise after rounding and nudge the diagonal so positive definiteness survives
        for i in [1..4] do for j in [i+1..4] do GI[j][i] := GI[i][j]; end for; end for;
        for i in [1..4] do GI[i][i] +:= 4; end for;
        L := LatticeWithGram(GI);
        SV := ShortVectors(L, Ceiling(SCL * Pbound * 1.000001));
        cfd := CC ! af(d);
        cf := (CC ! MoebiusMu(GCD(d, dat`N))) / (d * cfd);
        for sv in SV do
            v := Eltseq(sv[1]);
            al := &+[ v[i]*bs[i] : i in [1..4] ];
            t := Norm(al);
            if t le 0 then continue; end if;
            n := d * t;
            if not IsIntegral(n) then
                error "d*nrd(alpha) not integral: ", n;
            end if;
            n := Integers() ! n;
            Mx := &+[ (CC ! v[i]) * Ms[i] : i in [1..4] ];
            X := ((Mx[1][1] + Mx[2][2]) + Parent(CC.1).1 * (Mx[1][2] - Mx[2][1])) / 2;
            aX := Abs(X);
            if 4*pi*aX gt cut then continue; end if;
            w := WkFun(k, 4*pi*aX, RR);
            total +:= 2 * (CC ! af(n)) * cf * (X/aX)^k * w;   // 2: alpha and -alpha, k even
            nterms +:= 1;
        end for;
        if Verbose then printf "   d=%o: %o lattice points used\n", d, nterms; end if;
    end for;
    return total, nterms;
end function;
