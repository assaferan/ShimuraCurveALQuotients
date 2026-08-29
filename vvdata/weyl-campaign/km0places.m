// Where does the m=0 term at a nonzero isotropic coset actually differentiate?
//
// Schofer Thesis, Lemma 2.18 proof: for mu != 0, phi_mu(0) = 0, so
//     E_0(tau,s;phi_mu,+1) = -sqrt(pi) i v^(-s/2) Gamma((s+1)/2)/Gamma(s/2+1)
//                            * prod_p W_{0,p}(s,phi_mu),
// the archimedean factor is -i*pi at s=0 (nonzero), and the value vanishes by
// incoherence.  So kappa_mu(0) = -i*pi * (prod_p W_{0,p})'(0), and the log that
// appears is log p0 for whichever place p0 supplies the zero.
//
// km0poly.m showed the level prime is never that place: W_{0,N}(X,mu) = 1
// identically.  This script asks, on REAL bases at firing discriminants rather
// than a synthetic plane, which place does vanish -- and hence which logarithm
// the corrected m=0 term can produce.
//
// d/ds = -log(p) * X * d/dX, so at X=1 the s-derivative of a local factor is
// -log(p) * W'(1).  We report W'(1) and name the prime.
AttachSpec("ShimuraQuotients.spec");

for b in [<10,3,-8>, <15,2,-7>, <6,5,-19>] do
    D := b[1]; N := b[2]; d := b[3];
    Ld := ShimuraCurveLattice(D,N);
    Q := ChangeRing(Ld`Q, Integers()); Qrat := ChangeRing(Q, Rationals());
    lam := ElementOfNorm(Ld`Q, -d, Ld`O, Ld`basis_L);
    Lminus := Kernel(Transpose(Matrix(lam*Q)));
    Bm := ChangeRing(BasisMatrix(Lminus), Rationals());
    G := Bm*Qrat*Transpose(Bm);            // bilinear Gram of L_-
    Delta := Determinant(G);
    printf "\n===== X0^%o(%o), d = %o : disc(L_-) = %o = %o, |L_-^v/L_-| = %o =====\n",
           D, N, d, Delta, Factorization(Integers()!Delta), Abs(Integers()!Delta);
    printf "      (firing: N does not divide d.  N = %o)\n", N;

    // cosets of L_-^v / L_-  as  (a,b) * G^-1 * Bm,  (a,b) in Z^2 / G Z^2
    Gi := G^(-1);
    Gz := ChangeRing(G, Integers());
    Sm, P, R := SmithForm(Gz);             // Sm = P*Gz*R, invariant factors on diagonal
    e1 := Sm[1][1]; e2 := Sm[2][2];
    primes := [ p : p in PrimeDivisors(2*Abs(Integers()!Delta)) ];
    printf "      places examined: %o\n", primes;

    // NOTE.  At m=0 the local Whittaker is a POLYNOMIAL only where the defining
    // series terminates; where mu_p = 0 it need not, and this code path then
    // throws rather than returning a value.  We catch that and mark the place
    // NOT COMPUTED.  Do NOT read those as nonzero: at p | D the plane is
    // ANISOTROPIC, where Q(x)=0 has only the trivial solution, so such a factor
    // is a candidate to vanish and thereby to be the place of differentiation.
    // Deciding them needs the m=0 local density at a scaled anisotropic plane,
    // which is the next computation and is not done here.
    niso := 0; summary := [];
    for i in [0..Abs(e1)-1] do for j in [0..Abs(e2)-1] do
        if i eq 0 and j eq 0 then continue; end if;
        coef := Vector([Rationals()| i, j]) * ChangeRing(R^(-1), Rationals());
        cL := coef * Gi;                              // coordinates w.r.t. basis of L_-
        cL := Vector([Rationals()| x - Floor(x) : x in Eltseq(cL)]);   // reduce mod L_-
        if IsZero(cL) then continue; end if;
        mu := cL * Bm;                                // canonical small representative
        qmu := (mu*Qrat, mu)/2;
        if not IsIntegral(qmu) then continue; end if; // keep ISOTROPIC cosets only
        niso +:= 1;
        vals := []; vanish := [];
        for p in primes do
            ok := true; W := 0;
            try
                W := LocalWhittakerPolynomial(Rationals()!0, p, mu, Lminus, Q);
            catch e
                ok := false;
            end try;
            if ok then
                v1 := Evaluate(W, 1);
                Append(~vals, <p, Sprint(W), Sprint(v1), Sprint(Evaluate(Derivative(W),1))>);
                if v1 eq 0 then Append(~vanish, p); end if;
            else
                Append(~vals, <p, "non-terminating series (NOT COMPUTED)", "?", "?">);
            end if;
        end for;
        Append(~summary, vanish);
        if niso le 3 then
            printf "  mu = %o   Q(mu) = %o\n", Eltseq(mu), qmu;
            for t in vals do
                printf "      p=%-3o W(X) = %-38o W(1) = %-8o W'(1) = %-8o%o\n",
                       t[1], t[2], t[3], t[4],
                       t[1] eq N select "   <-- LEVEL PRIME" else
                        (t[3] eq "0" select "   <-- VANISHES" else "");
            end for;
            printf "      vanishing places: %o  ==> %o\n",
                   vanish, #vanish eq 1 select Sprintf("differentiates at log %o", vanish[1])
                           else (#vanish eq 0 select "no zero among the places we could compute -- the rest are UNDETERMINED"
                                 else "two or more zeros: derivative VANISHES");
        end if;
    end for; end for;
    printf "      isotropic nonzero cosets: %o\n", niso;
    printf "      vanishing-set multiset over all of them: %o\n", {* Sprint(v) : v in summary *};
end for;
printf "\nDONE\n"; quit;
