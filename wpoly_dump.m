// Dump the FULL local Whittaker POLYNOMIAL W_p(m, x) (x = p^{-s}) for every principal-part index of
// every Borcherds form of a base -- not just its value at x = 1, which is all LocalWhittakerAtOne
// returns.  The s-derivative is what produces the log p terms, so W_p'(1) is the quantity the m=0
// investigation needs and the pipeline currently throws away.
//
// usage: magma -b DD:=<D> NN:=<N> OUTNAME:=<file> wpoly_dump.m
// Emits, per (form, block, m, p):  [WP] the coefficient sequence of the polynomial, its value at 1,
// and its derivative at 1.

AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 0);

outname := "wpoly_out.txt";
if assigned OUTNAME then outname := OUTNAME; end if;
fh := Open(outname, "w");
emit := procedure(s) fprintf fh, "%o\n", s; Flush(fh); end procedure;

D := StringToInteger(DD); N := StringToInteger(NN);
Ld := ShimuraCurveLattice(D, N);
Q := Ld`Q; disc_grp := Ld`disc_grp; to_disc := Ld`to_disc; denom := Ld`denom;
M := IsOdd(D*N) select 4*D*N else 2*D*N;
Qint := ChangeRing(Q, Integers());
negQ := -Qint;
dQ := Determinant(Qint);
detprimes := Set(PrimeDivisors(dQ));
Lfull := RSpaceWithBasis(IdentityMatrix(Integers(), 3));

emit("START");
emit(Sprintf("=========== X0^%o(%o) ===========", D, N));
emit(Sprintf("det Q = %o = %o, M = %o", dQ, Factorization(dQ), M));

// the same bucketing as m0_multiplier / SchoferFormula0
mod_M_to_vecs := AssociativeArray([0..M-1]);
for j in [0..M-1] do mod_M_to_vecs[j] := []; end for;
i0 := 0;
for eta in disc_grp do
    if IsZero(eta) then i0 := eta; end if;
    v := ChangeRing(eta@@to_disc, Rationals());
    res := M*((v*ChangeRing(Q,Rationals()), v)/(2*denom^2));
    if not IsIntegral(res) then continue; end if;
    Append(~mod_M_to_vecs[Integers()!res mod M], eta);
end for;

dump_one := procedure(tag, k, eta, r, c)
    w_eta := ChangeRing(eta@@to_disc, Rationals())/denom;
    Sc := Sort([p : p in detprimes join Set(PrimeDivisors(Numerator(r)))]);
    for p in Sc do
        // unscaled polynomial, exactly as LocalWhittakerAtOne uses it before Evaluate(.,1)
        Wpol := WhittakerPolynomial(r, p, Vector(Rationals(), Eltseq(w_eta)), Lfull, negQ);
        R := Parent(Wpol);
        val := Evaluate(Wpol, 1);
        der := Evaluate(Derivative(Wpol), 1);
        emit(Sprintf("[WP] form=%o %o r=%o c=%o p=%o deg=%o coeffs=%o W(1)=%o Wprime(1)=%o",
                     k, tag, r, c, p, Degree(Wpol), Coefficients(Wpol), val, der));
    end for;
end procedure;

Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
ks := Sort([k : k in Keys(fs)]);
emit(Sprintf("keys = %o", ks));
for k in ks do
    foo := qExpansionAtoo(fs[k], 1);
    f0  := qExpansionAt0(fs[k], 1);
    for mm in [1..-Valuation(foo)] do
        c := Coefficient(foo, -mm);
        if c ne 0 then dump_one("oo", k, i0, Rationals()!mm, c); end if;
    end for;
    for j in [1..-Valuation(f0)] do
        c := Coefficient(f0, -j);
        if c eq 0 then continue; end if;
        r := (Rationals()!j)/M;
        for eta in mod_M_to_vecs[j mod M] do
            dump_one(Sprintf("zero j=%o eta=%o", j, Eltseq(ChangeRing(eta@@to_disc, Rationals())/denom)),
                     k, eta, r, c);
        end for;
    end for;
end for;
emit("DONE");
exit;
