// An exact model of a HYPERELLIPTIC star quotient X_0^D(N)/W, of any genus, by the recipe of
// Guo-Yang, Compositio 153 (2017), Example 43 (reproduced against the paper in
// tests/GuoYangExample43.m).  The genus-one case has its own driver, tools/genus1star.m.
//
//     magma -b D:=302 N:=1 tools/hyperellipticstar.m < /dev/null > hyp_302_1.log 2>&1
//
// Options (name:=value): stage (1 the search for the degree-two divisors, the default; 2 the model,
// which needs divisors), bd (class-number bound of the CM pool, default 2), divisors (the four
// divisors of stage 2, in Magma syntax: the first two are the fibres of the degree-two coordinate
// x, the last two those of the second coordinate y), verb (verbosity, default 1).
//
// STAGE 1 asks which degree-two divisors supported on the rational CM points are principal.  A
// function of degree two exists exactly when the curve is hyperelliptic, so a form found proves
// that (given that the form is a function -- Assumption 40, or its replacement when the Jacobian is
// torsion-free), and none found over an exhausted range is evidence the other way.  The divisors
// are enumerated in the general shape P_a + P_b - P_c - P_d, poles distinct or equal: a star whose
// hyperelliptic coordinate has two distinct poles is invisible to a search over double poles alone.
//
// STAGE 2 needs four forms, as in the paper: two whose divisors are fibres of x, so that the pair
// |x| and |x - t| fixes the sign of x at every CM point, and two more that do the same for y.  It
// prints the relation, the model, and the comparison that decides whether the model is right: the
// number of points of the curve over F_p against the trace formula, which needs no model.  THE FIT
// IS NOT SELF-CERTIFYING -- the values at the CM points do not determine the relation on their own
// and the divisor of y supplies the rest, so a wrong divisor yields a wrong curve of the right
// genus through the same points.  Read the trace-formula line before believing the model.
AttachSpec("ShimuraQuotients.spec"); SetColumns(0);
// ⚠ NMZ_TIMEOUT.  Normaliz is given 1800 s per polytope by default, which is right for the levels
// the hauptmodul search uses and wrong for these: three searches died on it in one day, at
// (4004, 1252), (8260, 2592) and (1540, 480, 343), each after hours of other work, and the rung
// they stopped on was reachable with more time.  A failure is loud and names the fix, so nothing
// is silently wrong, but the run is wasted; the warning is here because the setting belongs to
// the caller's environment and cannot be set from inside Magma.
if GetEnv("NMZ_TIMEOUT") eq "" then
    printf "⚠ NMZ_TIMEOUT is unset, so Normaliz gets 1800 s per polytope.  At these levels a deep\n"
         * "  rung needs hours: export NMZ_TIMEOUT=21600 before a long run.\n";
end if;
D := StringToInteger(D); N := StringToInteger(N);
stage := assigned stage select StringToInteger(stage) else 1;
bd := assigned bd select StringToInteger(bd) else 2;
verb := assigned verb select StringToInteger(verb) else 1;
SetVerbose("ShimuraQuotients", verb);

W := {m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1};
Xstar := CreateShimuraQuot(D, N, W);
Xstar`g := GenusShimuraCurveQuotient(D, N, W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
printf "X_0^%o(%o)*: genus %o\n", D, N, star`g;
error if star`g lt 2, "the star has genus < 2; for genus one use tools/genus1star.m";

// NONFUND admits the points by non-maximal orders, here and in the search's own lookup
nonfund := GetEnv("NONFUND") ne "";
rat, quad := RationalandQuadraticCMPoints(star : bd := bd, coprime_to_level := false,
                                                 AllOrders := nonfund);
ds := Reverse(Sort([p[1] : p in rat]));
printf "rational CM points (class number <= %o%o): %o\n", bd, nonfund select ", non-maximal orders included" else "", ds;
printf "quadratic CM cycles: %o\n", Sort([p[1] : p in quad]);

if stage eq 1 then
    // the effective divisors of degree two: two distinct rational points, one doubled rational
    // point, or one quadratic cycle
    eff := [[<ds[i], 1>, <ds[j], 1>] : j in [i+1..#ds], i in [1..#ds]]
       cat [[<d, 2>] : d in ds]
       cat [[<q[1], 1>] : q in quad];
    // a form with divisor A - B gives one with B - A by inversion, so only one of each pair is asked
    divisors := [];
    for i in [1..#eff] do
        for j in [1..#eff] do
            if j eq i then continue; end if;
            if not IsEmpty({c[1] : c in eff[i]} meet {c[1] : c in eff[j]}) then continue; end if;
            if j lt i then continue; end if;
            Append(~divisors, eff[i] cat [<c[1], -c[2]> : c in eff[j]]);
        end for;
    end for;
    printf "searching %o degree-two divisors over %o effective divisors\n", #divisors, #eff;
    // maxpole caps the ladder of pole orders at infinity: with it the search makes one pass at the
    // natural window and one at the given order, instead of climbing six rungs.  At a level with
    // many divisors the deep rungs are unaffordable -- enumeration costs about ten times as much
    // per further 22 in the rung -- so a capped pass is the bounded experiment: a form living
    // shallow is found cheaply, and nothing found means nothing shallow, not nothing.
    po := assigned maxpole select [StringToInteger(maxpole)] else [];
    t0 := Realtime();
    fs := BorcherdsForms(star, curves : Divisors := divisors, PoolClassNumber := bd,
                                        PoleOrders := po,
                                        IntegralSolution := true, AllowMissing := true);
    printf "done in %o s: %o of %o divisors carry a form\n", Realtime() - t0, #Keys(fs), #divisors;
    for k in Sort([k : k in Keys(fs)]) do
        printf "   divisor %o   pole order at oo %o\n", divisors[k], -Valuation(qExpansionAtoo(fs[k], 1));
    end for;
    exit;
end if;

error if not assigned divisors, "stage 2 needs divisors:=[...]";
divs := eval divisors;
error if #divs ne 4, "stage 2 needs four divisors: two fibres of x, then two of y";
fs := BorcherdsForms(star, curves : Divisors := divs, PoolClassNumber := bd,
                                    IntegralSolution := true);
for i in [1..4] do
    g := qExpansionAtoo(fs[i], 1);
    printf "form %o: divisor %o, pole order at oo %o, integral principal part %o\n",
           i, Set(DivisorOfBorcherdsForm(fs[i], star)), -Valuation(g),
           forall{m : m in [Valuation(g)..0] | IsIntegral(Coefficient(g, m))};
end for;

zeros := func<dv | Sort([c[1] : c in dv | c[2] gt 0])>;
vals := BorcherdsFormValues([fs[i] : i in [1..4]], ds, star);
xs := CoordinateFromBorcherdsForms(vals[1], vals[2], ds, zeros(divs[1]), zeros(divs[2]));
ys := CoordinateFromBorcherdsForms(vals[3], vals[4], ds, zeros(divs[3]), zeros(divs[4]));
printf "d  %o\nx  %o\ny  %o\n", ds, xs, ys;

F := HyperellipticRelation(ds, xs, ys, divs[3]);
printf "relation %o\n", F;
H, num, den := HyperellipticModelFromRelation(F);
printf "model %o\ngenus %o (the star has genus %o)\n", H, Genus(H), star`g;

// the comparison that decides it: the trace formula needs no model
f := HyperellipticPolynomials(IntegralModel(H));
agree := 0; tried := 0;
for p in [p : p in PrimesInInterval(3, 60) | (D*N) mod p ne 0] do
    if exists{c : c in Coefficients(f) | Valuation(Denominator(c), p) gt 0} then continue; end if;
    fp := PolynomialRing(GF(p))![GF(p) | c : c in Coefficients(f)];
    if Degree(fp) lt 2*star`g + 1 or not IsSquarefree(fp) then continue; end if;
    tried +:= 1;
    if #Points(HyperellipticCurve(fp)) eq ComputePointsViaTrace(Xstar, p, 1) then agree +:= 1; end if;
end for;
printf "trace formula: the model agrees at %o of %o good primes below 60\n", agree, tried;
printf "RESULT := <%o, %o, %o, %o>;\n", D, N, F, agree eq tried and tried ge 10;
