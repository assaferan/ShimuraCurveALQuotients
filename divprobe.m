// THEORY PROBE: is the m=0 multiplier sum_eta c_eta(0) recoverable from the Borcherds DIVISOR?
//
// Motivation: the Eisenstein/local-factor route to b_eta(m) is now known to be wrong in a way that is
// not a rescale (X0^10(11): old 4/5, 8/5, 12/5 -> true 1, 2, 1). But sum_eta c_eta(0) = 2 * weight(Psi_F),
// and the divisor of Psi_F is computed INDEPENDENTLY from the principal part (Guo-Yang Lemma 25, already
// implemented as DivisorOfBorcherdsForm). If Psi_F is a genuine weight-k form on X*, then
//     deg(div Psi) = k * deg(omega)
// so the divisor SOLVES for k without ever touching a local Whittaker factor. Not circular.
//
// This dumps, per form: the divisor <d, mult>, several candidate degree functionals, and (printed
// alongside by the caller) the measured multiplier, so the relation can be read off.
// usage: magma -b DD:=<D> NN:=<N> divprobe.m

AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 0);

D := StringToInteger(DD); N := StringToInteger(NN);
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);

ells := NumberOfEllipticPointsByCMOrder(star);
printf "[DIV] base %o %o  genus=%o  elliptic-orders=%o\n", D, N, star`g,
       Sort([q : q in Keys(ells)]);

for k in Sort([kk : kk in Keys(fs)]) do
    dv := DivisorOfBorcherdsForm(fs[k], star);
    if IsEmpty(dv) then
        printf "[DIV] %o %o form=%o EMPTY\n", D, N, k;
        continue;
    end if;
    smult := &+[pr[2] : pr in dv];
    // weight each disc by its CM-point count on X* (class number of the order, and the
    // number of optimal embeddings, are the natural candidates)
    wh := Rationals()!0; wn := Rationals()!0;
    for pr in dv do
        d := pr[1]; mu := pr[2];
        S := QuadraticOrder(BinaryQuadraticForms(d));
        h := ClassNumber(S);
        nemb := NumberOfOptimalEmbeddings(S, D, N);
        wh +:= mu * h;
        wn +:= mu * nemb;
    end for;
    printf "[DIV] %o %o form=%-4o div=%o  sum_mult=%o  sum_mult*h=%o  sum_mult*nemb=%o\n",
           D, N, k, dv, smult, wh, wn;
end for;
exit;
