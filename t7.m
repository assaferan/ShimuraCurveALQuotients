// The deliverable run.  For every Borcherds form of a base:
//   (A) check the numeric principal part of F_f against the pipeline's own two q-expansions
//       (GY Lemma 24 in the code's normalisation)  -- the correctness gate;
//   (B) report c_0(0) and sum over isotropic eta of c_eta(0).
// The slash points z_{t,gamma} and the eta values there do not depend on the form, so they are
// computed once and shared -- that is the difference between minutes and hours.
AttachSpec("ShimuraQuotients.spec");
Attach("vvlib.m");

D := 15; N := 2;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;
PREC := 120;  if assigned PR then PREC := StringToInteger(PR); end if;
KS := 64;     if assigned KK then KS := StringToInteger(KK); end if;
YY := 1;      if assigned YV then YY := StringToInteger(YV); end if;

t0 := Cputime();
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
ks := Sort([k : k in Keys(fs)]);
Rq := Parent(fs[ks[1]]); M := Rq`M; ds := Rq`ds;
printf "X0^%o(%o): %o forms, M = %o, ds = %o, keys = %o  (%o s)\n", D, N, #ks, M, ds, ks, Cputime(t0);

Ld := ShimuraCurveLattice(D, N);
CC := ComplexField(PREC);  ii := CC.1;  twopii := 2*Pi(CC)*ii;

t1 := Cputime();
S, Td, elts, i0 := VVRho(Ld, CC : Dual := true);
n := #elts;
for i in [1, 2, n] do
    assert Abs(&+[S[i][k]*ComplexConjugate(S[i][k]) : k in [1..n]] - 1) lt 1e-20;
end for;
printf "rho: |G| = %o  (%o s)\n", n, Cputime(t1);

Q := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
vs := [ChangeRing(g@@Ld`to_disc, Rationals()) : g in elts];
nmv := [ (vs[i]*Q, vs[i])/(2*dn^2) : i in [1..n] ];
res := [ IsIntegral(M*nmv[i]) select (Integers()!(M*nmv[i]) mod M) else -1 : i in [1..n] ];
iso := [i : i in [1..n] | IsIntegral(nmv[i])];
printf "#isotropic = %o ; #eta with M*nm not integral = %o\n", #iso, #[i : i in [1..n] | res[i] eq -1];

reps := VVCosetReps(M);
words := [VVSTWord(g) : g in reps];
assert &and[VVWordMatrix(words[k]) eq reps[k] : k in [1..#reps]];
t2 := Cputime();
U := ZeroMatrix(CC, #reps, n);
for k in [1..#reps] do
    v := VVRhoInvE0(S, Td, words[k], i0);
    for i in [1..n] do U[k][i] := v[i][1]; end for;
end for;
printf "%o coset rho-vectors (%o s)\n", #reps, Cputime(t2);

// ---- form-independent transcendental work: the slash point and eta values per (t, coset) ----
t3 := Cputime();
taus := [CC!(t/KS) + CC!YY*ii : t in [0..KS-1]];
FACT := ZeroMatrix(CC, KS, #reps);
ETAB := [[[CC | 0 : d in ds] : k in [1..#reps]] : t in [1..KS]];
maxIm := 0.0; minIm := 10.0^10;
for t in [1..KS] do
    for k in [1..#reps] do
        z := taus[t]; fac := CC!1;
        w := words[k];
        for r := #w to 1 by -1 do
            if w[r][1] eq "S" then fac /:= Sqrt(z); z := -1/z; else z := z + w[r][2]; end if;
        end for;
        FACT[t][k] := fac;
        ETAB[t][k] := [DedekindEta(d*z) : d in ds];
        maxIm := Maximum(maxIm, Im(z)); minIm := Minimum(minIm, Im(z));
    end for;
end for;
printf "eta table: %o points, Im in [%o, %o]  (%o s)\n", KS*#reps, minIm, maxIm, Cputime(t3);

slashval := function(f, t, k)
    tot := CC!0;
    for r in Exponents(f) do
        term := CC!(f`coeffs[r]);
        for i->d in ds do if r[i] ne 0 then term *:= ETAB[t][k][i]^r[i]; end if; end for;
        tot +:= term;
    end for;
    return FACT[t][k]*tot;
end function;

printf "\n%-6o %-7o %-7o %-13o %-24o %o\n", "form", "v_oo", "v_0", "max|PPerr|", "c_0(0)", "sum_iso c_eta(0)";
for kk in ks do
    f := fs[kk];
    foo := qExpansionAtoo(f, 80);  f0 := qExpansionAt0(f, 80);
    tt := Cputime();
    A := ZeroMatrix(CC, KS, #reps);
    for t in [1..KS] do for k in [1..#reps] do A[t][k] := slashval(f, t, k); end for; end for;
    Fv := A*U;                                  // Fv[t][i] = F_eta_i(tau_t)
    coef := function(i, nn)
        s := CC!0;
        for t in [1..KS] do s +:= Fv[t][i] * Exp(-twopii*CC!nn*taus[t]); end for;
        return s/KS;
    end function;
    // (A) principal part vs the code's prediction
    err := 0.0;  v0 := Minimum(0, Valuation(f0));
    for i in [1..n] do
        for j in [1..-v0] do
            if (j mod M) ne res[i] then continue; end if;
            nn := -(Rationals()!j)/M;
            pred := CC!Coefficient(f0, -j);
            if i eq i0 and IsIntegral(nn) then pred +:= Coefficient(foo, Integers()!nn); end if;
            err := Maximum(err, Abs(coef(i,nn) - pred));
        end for;
    end for;
    for nn in [Valuation(foo)..-1] do           // integer n beyond the 0-block's range
        if (-M*nn) le -v0 then continue; end if;
        err := Maximum(err, Abs(coef(i0,nn) - CC!Coefficient(foo, nn)));
    end for;
    c00 := coef(i0, 0);
    ciso := &+[CC | coef(i, 0) : i in iso];
    printf "%-6o %-7o %-7o %-13o %-24o %o   (%o s)\n", kk, Valuation(foo), Valuation(f0),
           ChangePrecision(err, 6), ChangePrecision(c00, 12), ChangePrecision(ciso, 12), Cputime(tt);
end for;
printf "\nTOTAL %o s\n", Cputime(t0);
quit;
