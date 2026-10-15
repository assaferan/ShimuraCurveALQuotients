// Diagnostic for the Atkin-Lehner layer at D = 6, N = 1.  dim M_4(Gamma) = 1, and the Atkin-Lehner
// group acts on that line, so each w_m must act by +1 or -1.  Solve the Gamma-only system for the
// weight-4 form, then measure the residual of the Atkin-Lehner relation and the implied eigenvalue.
// If the relation is implemented correctly the residual is at truncation level for one sign.
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vw3.m";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vwal.m";

PREC := 100; NT := 50; Q := 300; CC := ComplexField(PREC);
D := 6; N := 1;
A := QuaternionAlgebra(D); O := QuaternionOrder(A, N);
G := FuchsianGroup(O); DD := UnitDisc( : Precision := PREC);
FD := FundamentalDomain(G, DD);
rho := Max([Abs(ComplexValue(v)) : v in FD]);
H := UpperHalfPlane();
p := ComplexValue(DiscToPlane(H, DD ! 0) : Precision := PREC);
mrep := FuchsianMatrixRepresentation(A : Precision := PREC);
sampG := ReducedSamples(G, O, DD, mrep, p, rho, Q, CC);
printf "rho %o, %o Gamma samples, NT %o, max |Gamma target| %o\n", RealField(6)!rho, #sampG, NT, RealField(4)!Max([Abs(s[2]) : s in sampG]);

k := 4;
MG := NormaliseColumns(AutSystem(sampG, k, NT, CC));
svG := SingVals(MG);
printf "weight %o Gamma-only: sv[1..4] = %o\n", k, [RealField(3)!s : s in svG[1..4]];
vec, s0 := SmallestVec(MG);
// undo the column scaling to recover the actual coefficients b_n
Mraw := AutSystem(sampG, k, NT, CC);
nrm := [ Sqrt(&+[ Abs(Mraw[i][j])^2 : i in [1..Nrows(Mraw)] ]) : j in [1..NT] ];
b := [ nrm[j] gt 0 select vec[j]/nrm[j] else vec[j] : j in [1..NT] ];
bmax := Max([Abs(x) : x in b]);
printf "coefficient vector: max |b_n| = %o, smallest sv = %o\n", RealField(3)!bmax, RealField(3)!s0;

// residual of the Gamma relation on the recovered form, as the scale to beat
function Resid(samples, b, k)
    worst := 0.0;
    for s in samples do
        w := s[1]; wp := s[2]; jk := DiscJ(s[3], w)^k;
        lhs := &+[ b[n+1]*wp^n : n in [0..#b-1] ];
        rhs := jk * &+[ b[n+1]*w^n : n in [0..#b-1] ];
        sc := Max(Abs(lhs), Abs(rhs));
        if sc gt 0 then worst := Max(worst, Abs(lhs-rhs)/sc); end if;
    end for;
    return worst;
end function;
printf "Gamma relation residual on the recovered form: %o\n", RealField(3)!Resid(sampG, b, k);

for m in [2, 3, 6] do
    ok, mu := ALElement(O, m, 3);
    printf "--- m = %o, mu = %o, nrd = %o ---\n", m, mu, Norm(mu);
    sampA := ALSamples(G, O, DD, mrep, p, rho, Q, CC, [mu]);
    printf "   %o Atkin-Lehner samples\n", #sampA;
    // consistency of the composite: target must equal the matrix applied to the source
    worstc := 0.0;
    for s in sampA do worstc := Max(worstc, Abs(s[2] - ActOn(s[3], s[1]))); end for;
    printf "   composite consistency (target vs matrix.source): %o\n", RealField(3)!worstc;
    // modulus range actually sampled
    printf "   |source| %o, max |target| %o\n", RealField(4)!Abs(sampA[1][1]), RealField(4)!Max([Abs(s[2]) : s in sampA]);
    // is the composite in SU(1,1)?  |a|^2-|c|^2 should be 1 and b = conj(c), d = conj(a)
    Wt := sampA[1][3];
    printf "   SU(1,1) check: |a|^2-|c|^2 = %o, |b-conj(c)| = %o, |d-conj(a)| = %o\n",
        RealField(3)!(Abs(Wt[1][1])^2 - Abs(Wt[2][1])^2),
        RealField(3)!Abs(Wt[1][2] - ComplexConjugate(Wt[2][1])),
        RealField(3)!Abs(Wt[2][2] - ComplexConjugate(Wt[1][1]));
    // eigenvalue: compare both signs of the relation
    // how far the reduction actually got: modulus before and after
    Wmu := DiscMatrix(mrep, mu, CC ! p, CC);
    pre := [ Abs(ActOn(Wmu, s[1])) : s in sampA ];
    printf "   |mu.w| before reduction: max %o min %o ; |target| after: max %o, over rho: %o of %o\n",
        RealField(4)!Max(pre), RealField(4)!Min(pre),
        RealField(4)!Max([Abs(s[2]) : s in sampA]),
        #[s : s in sampA | Abs(s[2]) gt rho + 1e-8], #sampA;
    for sgn in [1, -1] do
        worst := 0.0;
        for s in sampA do
            if Abs(s[2]) gt rho + 1e-8 then continue; end if;
            w := s[1]; wp := s[2]; jk := DiscJ(s[3], w)^k;
            lhs := &+[ b[n+1]*wp^n : n in [0..#b-1] ];
            rhs := sgn * jk * &+[ b[n+1]*w^n : n in [0..#b-1] ];
            sc := Max(Abs(lhs), Abs(rhs));
            if sc gt 0 then worst := Max(worst, Abs(lhs-rhs)/sc); end if;
        end for;
        printf "   residual (well-reduced samples only) with sign %o : %o\n", sgn, RealField(3)!worst;
    end for;
end for;
quit;
