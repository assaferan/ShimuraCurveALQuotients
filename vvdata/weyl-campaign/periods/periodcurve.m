// A GENUS-TWO CURVE FROM A MODULAR ABELIAN SURFACE, BY ITS PERIODS.
//
//   magma -b LEVEL:=23 DIGITS:=200 periodcurve.m < /dev/null
//   magma -b LEVEL:=2002 DIGITS:=300 WANT:=13 periodcurve.m < /dev/null   (WANT: pick the factor
//                                                                          with this RM field disc)
//
// WHY THIS ROUTE.  The Jacobian of a genus-two Atkin-Lehner star is a modular abelian surface, and
// its periods come from the classical newform, where precision is CHEAP: measured below, 198 real
// digits in well under a second.  The Nelson-Shimizu route on the Shimura curve costs minutes per
// evaluation and its recognition stalls near twelve digits, which is why the star models are stuck.
//
// THE STEP THAT IS THE ACTUAL CONTENT: THE POLARIZATION.  Magma's Periods and its
// IntersectionPairing are NOT in the same homology basis -- symplectifying the pairing and applying
// it to the periods leaves the Riemann matrix asymmetric by 0.94 -- so they cannot be composed.
// The polarization is therefore recovered FROM the periods: for a period matrix Pi (g x 2g, columns
// the lattice generators) and an integral alternating Q, Riemann's first relation is Pi Q Pi^t = 0,
// which for g = 2 is one complex equation in six unknowns.  Its integral solutions are the
// Neron-Severi group (rank 2 for a surface with real multiplication), found here by LLL; the
// principal polarizations among them are those of Pfaffian +-1, cut down by Riemann's second
// relation (positivity).  This is the computational face of Gonzalez-Guardia-Rotger, Acta Arith.
// 116 (2005).
//
// ⚠ PRECISION IS NOT WHAT Periods CLAIMS.  Its second argument is the NUMBER OF TERMS, not digits:
// measured at level 23, accuracy runs at about n/7 digits and saturates against the ambient complex
// field, which defaults to 30.  SetDefaultRealFieldPrecision lifts that, giving p - 2 real digits.
// Keying LLL's multiplier or any tolerance to the REQUESTED precision instead of the achieved one
// makes LLL fit noise and report six spurious "relations"; every threshold below is keyed to the
// accuracy actually measured, by the Riemann residual itself.
//
// CERTIFICATION, and why it needs no prior knowledge of the answer.  The curve is recognised from
// its invariants, so it must be tied back to the right newform: #C(F_p) = p + 1 - Tr(a_p) at many
// primes, with the quadratic twist chosen to fit.  That check is available at a level where the
// answer is unknown, which is the point.
//
// CONTROLS.  LEVEL 23: J_0(23) IS the Jacobian of the genus-two X_0(23), so the route must return
// it.  A second control at a level where the 2-dimensional factor is a PROPER quotient of J_0(N)
// is essential, because at 23 the polarization is principal for free.

AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
N := StringToInteger(LEVEL);
P := StringToInteger(DIGITS);
want := (assigned WANT) select StringToInteger(WANT) else 0;
R8 := RealField(8);
SetDefaultRealFieldPrecision(P);
g := 2;

A := JZero(N);
dec := Decomposition(A);
printf "=== J_0(%o): dimension %o, decomposition %o\n", N, Dimension(A), [ Dimension(B) : B in dec ];
twos := [ B : B in dec | Dimension(B) eq 2 ];
error if IsEmpty(twos), "no two-dimensional factor at this level";
// when several, pick by the discriminant of the real multiplication field
Af := twos[1];
if want ne 0 then
    hit := false;
    for B in twos do
        f := CharacteristicPolynomial(HeckeOperator(B, 3));
        d := 0;
        try d := SquarefreeFactorization(Discriminant(MaximalOrder(NumberField(f)))); catch e ; end try;
        if d eq want then Af := B; hit := true; break; end if;
    end for;
    error if not hit, Sprintf("no 2-dimensional factor with RM field of discriminant %o", want);
end if;
printf "factor: %o\n", Af;
printf "proper quotient of J_0(%o): %o\n", N, Dimension(A) ne 2;

terms := 8*P;
t0 := Cputime();
Pi := Transpose(Matrix([ Eltseq(v) : v in Periods(Af, terms) ]));
CC := BaseRing(Pi);
printf "periods: %o terms, ambient precision %o, %o s\n", terms, Precision(CC), R8 ! Cputime(t0);
error if Ncols(Pi) ne 2*g, "expected 2g period vectors";

pairs := [ <i,j> : j in [i+1..2*g], i in [1..2*g] ];
coef := [ Pi[1][p[1]]*Pi[2][p[2]] - Pi[1][p[2]]*Pi[2][p[1]] : p in pairs ];
scale := Maximum([ Abs(c) : c in coef ]);

// integer relations.  A genuine relation has residual at the periods' accuracy; a spurious one has
// residual of the order of its own coefficients, which is what chance gives -- so the test is the
// RATIO, not the residual.
// ⚠ THE MULTIPLIER MUST BE MODERATE, NOT MAXIMAL.  With M at the full precision the lattice is so
// skewed that LLL returns only its long vectors -- heads of 1e47 and above -- and the genuine
// short relations never surface: at level 39 and 67 that produced "no integral polarization"
// while principal ones exist.  So M is swept, and only relations that recur across several
// multipliers are kept, which is also the test that they are not artefacts of one scaling.
seen := AssociativeArray();
for e in [20, 30, 40, 60] do
    M := Round(10^e / scale);
    L := ZeroMatrix(Integers(), #pairs, #pairs + 2);
    for k in [1..#pairs] do
        L[k][k] := 1;
        L[k][#pairs+1] := Round(M * Real(coef[k]));
        L[k][#pairs+2] := Round(M * Imaginary(coef[k]));
    end for;
    Bm := LLL(L);
    for r in [1..Nrows(Bm)] do
        hd := [ Bm[r][k] : k in [1..#pairs] ];
        h := Maximum([ Abs(x) : x in hd ]);
        if h eq 0 or h gt 10^6 then continue; end if;
        tail := Maximum(Abs(Bm[r][#pairs+1]), Abs(Bm[r][#pairs+2]));
        if tail gt 1000*h then continue; end if;
        key := Sprint(hd lt [ -x : x in hd ] select hd else [ -x : x in hd ]);
        if IsDefined(seen, key) then seen[key] := seen[key] + 1; else seen[key] := 1; end if;
    end for;
end for;
cands := [];
for k in Keys(seen) do
    if seen[k] ge 2 then
        Append(~cands, eval k);
        printf "  relation %o  (recurs at %o of 4 multipliers)\n", eval k, seen[k];
    end if;
end for;
printf "%o stable integral relations; the Neron-Severi rank of an RM surface is 2\n", #cands;
error if IsEmpty(cands), "no integral polarization found";

function AltFromVec(v)
    Q := ZeroMatrix(Integers(), 2*g, 2*g);
    for k in [1..#pairs] do
        Q[pairs[k][1]][pairs[k][2]] := v[k]; Q[pairs[k][2]][pairs[k][1]] := -v[k];
    end for;
    return Q;
end function;
pf := func< Q | Q[1][2]*Q[3][4] - Q[1][3]*Q[2][4] + Q[1][4]*Q[2][3] >;

princ := [];
for cs in CartesianPower([-6..6], #cands) do
    v := [ &+[ cs[t]*cands[t][k] : t in [1..#cands] ] : k in [1..#pairs] ];
    if forall{ x : x in v | x eq 0 } then continue; end if;
    Q := AltFromVec(v);
    if Abs(pf(Q)) ne 1 then continue; end if;
    // NB a nested sequence, filled row by row: Magma's [ f(i,j) : i in I, j in J ] varies i
    // FASTEST, which fills a g x 2g matrix in the wrong order and silently scrambles the conjugate.
    cPi := Matrix(CC, g, 2*g, [[ ComplexConjugate(Pi[i][j]) : j in [1..2*g] ] : i in [1..g] ]);
    H := CC.1 * Pi * ChangeRing(Q, CC) * Transpose(cPi);
    herm := Matrix(RealField(P), g, g, [ Real((H[i][j]+H[j][i])/2) : i, j in [1..g] ]);
    if IsPositiveDefinite(herm) or IsPositiveDefinite(-herm) then Append(~princ, Q); end if;
end for;
printf "%o principal polarizations among small combinations\n", #princ;
error if IsEmpty(princ), "no principal polarization: the surface may admit none";

for Q in princ do
    // the Riemann residual at the accepted polarization measures the periods' REAL accuracy
    v := [ Q[p[1]][p[2]] : p in pairs ];
    res := Abs(&+[ v[k]*coef[k] : k in [1..#pairs] ]) / scale;
    acc := Round(-Log(10, res + 10^(-2*P)));
    printf "\npolarization %o : Riemann residual %o, so about %o real digits\n", v, R8 ! res, acc;

    E := Matrix(Integers(), 2*g, 2*g, [ Integers() ! x : x in Eltseq(Q^-1) ]);
    F, T := FrobeniusFormAlternating(E);
    Pin := Pi * Transpose(ChangeRing(T, CC));
    P1 := Submatrix(Pin,1,1,g,g); P2 := Submatrix(Pin,1,g+1,g,g);
    if Abs(Determinant(P1)) lt 10^(-acc div 2) then printf "  singular half, skipped\n"; continue; end if;
    tau := P1^-1 * P2;
    sym := Maximum([ Abs(tau[i][j]-tau[j][i]) : i, j in [1..g] ]);
    if sym gt 10^(-acc div 2) then printf "  Riemann matrix not symmetric (%o), skipped\n", R8 ! sym; continue; end if;
    tau := Matrix(CC, g, g, [ (tau[i][j]+tau[j][i])/2 : i, j in [1..g] ]);
    if not IsPositiveDefinite(Matrix(RealField(P), g, g, [ Imaginary(tau[i][j]) : i, j in [1..g] ])) then
        tau := -tau;
    end if;

    ok := true; rh := [];
    try rh := RosenhainInvariants(tau); catch e ok := false; printf "  Rosenhain refused: %o\n", e`Object; end try;
    if not ok then continue; end if;
    PC<x> := PolynomialRing(CC);
    f6 := x*(x-1)*&*[ x - e : e in rh ];
    // The Cardona-Quer absolute invariants, with Magma's own reconstruction from them.  Feeding
    // Igusa-Clebsch RATIOS to Mestre instead -- [1, j2/j1, j3/j1, 1/j1], denominators about 1e22
    // at level 2002 -- produces a model so large that reducing it ran 46 minutes and 19 GB.
    j := G2Invariants(HyperellipticCurve(f6));   // G2Invariants needs a curve, not a bare polynomial
    printf "  absolute invariants (real parts) %o\n", [ R8 ! Real(a) : a in j ];
    printf "  imaginary parts relative %o\n", [ R8 ! (Abs(Imaginary(a))/(1+Abs(a))) : a in j ];

    // recognition, with a margin: a rational of denominator q fits ANY real to about 1/q^2, so the
    // fit must beat that rate by a wide factor before it counts
    rec := []; allrec := true;
    for a in j do
        re := Real(a); bq := 0; bm := 0;
        for e in [1..(acc div 3)] do
            q := BestApproximation(re, 10^e);
            er := Maximum(Abs(RealField(P) ! q - re), 10^(-acc));
            ch := 1/(RealField(P) ! Maximum(Denominator(q),1))^2;
            if ch/er gt bm then bm := ch/er; bq := q; end if;
        end for;
        Append(~rec, bq);
        printf "    %o   margin %o\n", bq, R8 ! bm;
        if bm lt 10^6 then allrec := false; end if;
    end for;
    if not allrec then printf "  NOT RECOGNISED with a safe margin\n"; continue; end if;

    C := HyperellipticCurveFromG2Invariants([ Rationals() | r : r in rec ]);
    // Mestre's output is far from minimal -- 128-digit coefficients at level 23 -- and reducing a
    // non-minimal model mod small p gives bad or singular reductions, so the point counts below
    // would be meaningless without this step.
    try C := ReducedMinimalWeierstrassModel(IntegralModel(C)); catch e
        printf "  could not reduce the model: %o\n", e`Object; end try;
    printf "  CANDIDATE CURVE: %o\n", C;

    // CERTIFICATION.  Mestre's construction gives the curve only up to quadratic twist, and the
    // right twist is not generally small, so it is DETERMINED rather than guessed: at each good
    // prime the count is p + 1 -+ a_p, and the sign pattern is the quadratic character of the
    // twisting discriminant.  Recovering d from that pattern and then checking every prime is a
    // genuine test -- a wrong curve has no single d fitting all the signs.
    aps := AssociativeArray();
    for p in PrimesInInterval(3, 200) do
        // T_p acts on the 4-dimensional homology as a_p on K^2, so its trace is TWICE the trace
        // of a_p from K to Q, which is what p + 1 - #C(F_p) equals.  Verified against X_0(23):
        // the ratio is exactly 2 at every prime tested.
        if N mod p ne 0 then
            aps[p] := (Integers() ! Trace(HeckeOperator(Af, p))) div 2;
        end if;
    end for;
    ps := Sort([ k : k in Keys(aps) ]);
    sgn := AssociativeArray(); usable := [];
    for p in ps do
        ok2 := true; n := 0;
        try n := #Points(ChangeRing(C, GF(p))); catch e ok2 := false; end try;
        if not ok2 then continue; end if;
        // a_p = 0 makes the two signs indistinguishable, so such a prime carries no information
        // about the character and must not be fed to the fit -- it is still checked afterwards
        if aps[p] eq 0 then
            if n eq p + 1 then Append(~usable, p); end if;
            continue;
        end if;
        if n eq p + 1 - aps[p] then sgn[p] := 1; Append(~usable, p);
        elif n eq p + 1 + aps[p] then sgn[p] := -1; Append(~usable, p);
        end if;
    end for;
    printf "  point counts of the shape p+1-+a_p at %o of %o primes\n", #usable, #ps;
    if #usable lt 10 then
        printf "  TOO FEW: the curve does not have this newform's counts at all\n"; continue;
    end if;
    // The twisting discriminant can only ramify where the model has bad reduction, so d runs over
    // signed squarefree products of those primes -- a small set -- rather than over a range.
    bad := PrimeDivisors(Integers() ! Discriminant(C));
    printf "  bad primes of the model: %o\n", bad;
    printf "  sign pattern: %o\n", [ <q, sgn[q]> : q in usable | IsDefined(sgn, q) ];
    dfound := 0;
    for sub in Subsets({ x : x in bad }) do
        for sg in [1, -1] do
            d0 := sg * (IsEmpty(sub) select 1 else &*[ x : x in sub ]);
            dd := (d0 mod 4 in {0,1}) select d0 else 4*d0;
            if dd eq 0 then continue; end if;
            if forall{ q : q in usable | not IsDefined(sgn, q) or KroneckerSymbol(dd, q) eq sgn[q] }
            then dfound := dd; break; end if;
        end for;
        if dfound ne 0 then break; end if;
    end for;
    if dfound eq 0 then
        printf "  no single quadratic character fits the sign pattern: the curve is WRONG\n";
    else
        D := (dfound eq 1) select C else QuadraticTwist(C, dfound);
        allp := true; np := 0;
        for p in usable do
            if #Points(ChangeRing(D, GF(p))) ne p + 1 - aps[p] then allp := false; break; end if;
            np +:= 1;
        end for;
        printf "  twist by %o : %o\n", dfound,
               allp select Sprintf("CERTIFIED at all %o primes", np) else "FAILS after twisting";
        if allp then
            try printf "  model: %o\n", ReducedMinimalWeierstrassModel(IntegralModel(D));
            catch e printf "  (model resists reduction; the invariants are the certificate)\n"; end try;
        end if;
    end if;
end for;
quit;
