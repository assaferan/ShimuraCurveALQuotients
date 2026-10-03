// A conductor prime OUTSIDE the level: q = 7 on X_0^15(2) at d = -588 (= 49 * -12, 7 split in Q(sqrt -3),
// also conductor 2), -735 (= 49 * -15, 7 inert), -1960 (= 49 * -40, 7 split).  Table 45 is silent there,
// so the test is internal: the nine Borcherds forms are functions of ONE Hauptmodul s with known
// divisors (zeros at s = 0, 2, -1/12, 5/4), so with A = N|s|, B = N|s-2|, C = N|s+1/12|, E = N|s-5/4|
// (norms over the star points of discriminant d) and w_k = log|value_k| - npts log C_k:
//   w11 = w-1 + w-2,  w13 = w-1,  w15 = w-2,  w9 = w14 + w-2,  w10 = w14 + w-1,  w12 = w14 + w-1 + w-2.
// Checked WITH the q-term (the code) and WITHOUT it (subtracting M0FibreCorrection's value).
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
SetVerbose("ShimuraQuotients", 1);
D := 15; N := 2;
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
Ld := ShimuraCurveLattice(D, N);
ks := [-2, -1, 9, 10, 11, 12, 13, 14, 15];
etas := [fs[k] : k in ks];
Cs := AssociativeArray();   // the divisor constants, tests/M0PoleSum.m (3b), read at d = -7
Cs[-2] := 16/3; Cs[-1] := 320/3; Cs[9] := 20/3; Cs[10] := 100/3; Cs[11] := 1280/9; Cs[12] := 1600/9; Cs[13] := 320/3; Cs[14] := 5/4; Cs[15] := 16/3;
logC := func< c | &+[LogSum() ] + &+([LogSum()] cat [LogSum(Rationals()!t[2], t[1]) : t in Factorization(Numerator(c))] cat [LogSum(Rationals()!-t[2], t[1]) : t in Factorization(Denominator(c))]) >;
M := 2*D*N;
foos := [qExpansionAtoo(e, 1) : e in etas]; f0s := [qExpansionAt0(e, 1) : e in etas];
for d in [-588, -735, -1960] do
    d0 := FundamentalDiscriminant(d); _, f := IsSquare(d div d0);
    OK := MaximalOrder(QuadraticField(d)); O := sub<OK | f>;
    nd := NumberOfOptimalEmbeddings(O, D, N); npts := nd div 8;
    printf "\n=== d = %o = %o^2 * %o, h(R_f) = %o, optimal embeddings %o, star points %o\n", d, f, d0, PicardNumber(O), nd, npts;
    vals := SchoferFormula(etas, d, D, N, Ld : PointDegree := npts);
    lam := ElementOfNorm(Ld`Q, -d, Ld`O, Ld`basis_L);
    qcorr := M0FibreCorrection(foos, f0s, [Rationals() | 0 : e in etas], d, 7, Ld`Q : Lambda := lam, M := M, Unimodular := true);
    printf "q = 7 corrections per point (units of log 7): %o\n", qcorr;
    W := AssociativeArray(); W0 := AssociativeArray();
    for i->k in ks do
        W[k] := vals[i] - npts*logC(Cs[k]);
        W0[k] := W[k] - npts*qcorr[i]*LogSum(Rationals()!1, 7);
        printf "  form %-3o value %-30o\n", k, vals[i];
    end for;
    rels := [* <"w11 = w-1 + w-2", 11, [-1, -2]>, <"w13 = w-1", 13, [-1]>, <"w15 = w-2", 15, [-2]>,
              <"w9 = w14 + w-2", 9, [14, -2]>, <"w10 = w14 + w-1", 10, [14, -1]>, <"w12 = w14 + w-1 + w-2", 12, [14, -1, -2]> *];
    for tag in ["WITH q-term", "WITHOUT q-term"] do
        T := tag eq "WITH q-term" select W else W0;
        for r in rels do
            lhs := T[r[2]]; rhs := &+[T[k] : k in r[3]];
            printf "  %-15o %-24o lhs - rhs = %o\n", tag, r[1], lhs - rhs;
        end for;
    end for;
end for;
