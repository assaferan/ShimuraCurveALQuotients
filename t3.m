// Step 3: pin the metaplectic slash against the pipeline's VALIDATED cusp expansions.
//  (a) invariance:  |f|gamma| = |f| for gamma in Gamma_0(M)   (weight 1/2, unitary character)
//  (b) the S-slash:  compare (f|S~)(tau) against the code's qExpansionAt0 series, which the m>0
//      machinery is validated on.  Any constant ratio = the normalisation tying the two conventions.
AttachSpec("ShimuraQuotients.spec");
Attach("vvlib.m");
SetVerbose("User1", 0);

D := 15; N := 2;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;

t0 := Cputime();
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
ks := Sort([k : k in Keys(fs)]);
printf "built %o forms in %o s; keys = %o\n", #ks, Cputime(t0), ks;

f := fs[ks[1]];
R := Parent(f);
printf "R`M = %o   R`disc = %o   ds = %o\n", R`M, R`disc, R`ds;
M := R`M;

CC := ComplexField(80);
ii := CC.1;

// ---- (a) Gamma_0(M) invariance of |f| ----
Sm := Matrix(Integers(), 2, 2, [0,-1,1,0]);
gam := Matrix(Integers(), 2, 2, [1, 0, M, 1]);          // in Gamma_0(M)
gam2 := Matrix(Integers(), 2, 2, [1+M, 1, M, 1]);        // det = 1+M-M = 1, in Gamma_0(M)
for g in [gam, gam2] do
    w := VVSTWord(g);
    assert VVWordMatrix(w) eq g;
    for tau in [CC | 0.11 + 0.9*ii, -0.37 + 1.4*ii] do
        a := VVSlashEval(f, w, tau);
        b := VVEtaEval(f, tau);
        printf "(a) gamma=%o tau=%o   |f|gamma| / |f| = %o\n", Eltseq(g), tau, Abs(a)/Abs(b);
    end for;
end for;

// ---- (b) the S-slash vs qExpansionAt0 ----
ser0 := qExpansionAt0(f, 60);
printf "(b) val(qExpansionAt0) = %o   val(qExpansionAtoo) = %o\n",
       Valuation(ser0), Valuation(qExpansionAtoo(f,60));
// PREDICTION, from unwinding SAction against eta(-1/z) = sqrt(z/i) eta(z) termwise:
//   SAction divides the r-th monomial by scale_r = sqrt(prod d^{r_d} * disc) and multiplies by M;
//   the prod d^{r_d} cancels the M^k/sqrt(prod d^{r_d}) of the eta transformation, leaving
//        (f|S~)(M z) = e(-1/8) * sqrt(disc) / M  *  ser0(q),   q = e(z).
pred := Exp(-Pi(CC)*ii/4) * Sqrt(CC!R`disc) / M;
printf "(b) predicted constant e(-1/8)*sqrt(disc)/M = %o\n", pred;
wS := [<"S",0>];
for tau in [CC | 90*ii, 12.0 + 75*ii, -20.0 + 140*ii] do
    lhs := VVSlashEval(f, wS, tau);                       // (f|S~)(tau)
    q := Exp(2*Pi(CC)*ii*tau/M);                          // series variable is q^{1/M}
    rhs := &+[CC | Coefficient(ser0, n)*q^n : n in [Valuation(ser0)..55]];
    printf "(b) Im(tau)=%-6o  ratio = %o   ratio/pred = %o\n", Im(tau), lhs/rhs, (lhs/rhs)/pred;
end for;
quit;
