// Measure the outer-m=0 multipliers of the two Hauptmodul forms on a base, WITHOUT the numeric
// vector-valued oracle: drive the pipeline as far as the Hauptmodul table and hand it to
// HauptmodulM0Residuals (SchoferFormula.m), where the argument is written out.
//
//   magma DD:=14 NN:=5 vvdata/gtsweep.m
//
// Cost is seconds per base against ~20 minutes for the oracle, and it reaches bases the oracle
// cannot -- X0^10(11) has |disc grp| = 24200, where the oracle's 5-hour run failed its own gate.
// The limitation is that it sees rows -1 and -2 only, not the whole form list.
//
// Validated against every independently known answer: X0^6(5) -> (0,0) and X0^10(3) -> (0,0), both
// matching the oracle; X0^15(2) -> (4,2), matching the ground truth measured there; X0^14(5) -> (1,1),
// which the oracle then CONFIRMED to eleven digits (vvdata/run_14_5.log).
AttachSpec("ShimuraQuotients.spec");

D := 14; N := 5;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;
PRECF := 100; if assigned PF then PRECF := StringToInteger(PF); end if;

t0 := Cputime();
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};

// --- ValuesAtCMPoints(Xstar, curves) up to, but NOT including, find_signs_hauptmodul -------------
fs := BorcherdsForms(star, curves : Prec := PRECF);
d_divs := &cat[[T[1] : T in DivisorOfBorcherdsForm(f, star)] : f in [fs[-1], fs[-2]]];
must_use := Set(d_divs);

// The measurement needs discriminants whose firing status differs from the normalising points', and
// CandidateDiscriminants' coprime-to-level filter throws out exactly the N-divisible ones.  On
// X0^14(5) two survive as Hauptmodul anchors; elsewhere admit a few through the documented Keep hook
// (the mechanism tests/ExternalCMValues.m uses to reach a level-divisible discriminant).
rat_all := RationalandQuadraticCMPoints(star : coprime_to_level := false);
nonfiring := [ p[1] : p in rat_all
               | IsEmpty(PrimeDivisors(N div GCD(N, FundamentalDiscriminant(p[1])))) ];
extra := [ d : d in nonfiring | d notin must_use ];
extra := extra[1 .. Minimum(3, #extra)];
printf "X0^%o(%o): non-firing rational candidates %o; admitting %o beyond the anchors\n",
       D, N, nonfiring, extra;
must_use := must_use join Set(extra);

all_cm_pts := CandidateDiscriminants(star, curves : Keep := must_use);
abs_tab, all_cm_pts := AbsoluteValuesAtCMPoints(star, curves, all_cm_pts, fs
                                                : MaxNum := Maximum(7, #must_use + 2),
                                                  Prec := PRECF, Include := must_use);
ReduceTable(abs_tab);

ds    := abs_tab`Discs;
s     := [* RationalNumber(x) : x in abs_tab`Values[abs_tab`sIndex] *];
stild := [* RationalNumber(x) : x in abs_tab`Values[abs_tab`sTildeIndex] *];
degs  := [ Degree(abs_tab`FldsOfDefn[abs_tab`Xstar`CurveID][d][1]) : d in ds ];
printf "table built (%o s)\n  ds     = %o\n  s      = %o\n  stilde = %o\n  degs   = %o\n",
       Cputime(t0), ds, s, stild, degs;
for i->d in ds do
    printf "    d = %-8o deg %o  fires %o\n", d, degs[i],
           not IsEmpty(PrimeDivisors(N div GCD(N, FundamentalDiscriminant(d))));
end for;

sols := HauptmodulM0Residuals(s, stild, ds, degs, N);
if IsEmpty(sols) then
    printf "\nNOT MEASURABLE on this base: either no discriminant differs in firing status from the\n"
         * "normalising ones, or no multiple of log N explains the table.  Admit more N-divisible\n"
         * "discriminants through Keep and retry.\n";
    quit;
end if;
printf "\nadmissible residual pairs <r(-1), r(-2)>: %o\n", sols;
if #sols gt 1 then
    printf "NOT DETERMINED: %o pairs survive; more informative points are needed.\n", #sols;
    quit;
end if;

sol := Rep(sols);
res := AssociativeArray();  res[-1] := sol[1];  res[-2] := sol[2];
// The residual is r = true - applied.  TWO CORRECTIONS, 2026-08-29:
//
//  (1) `applied := IsEven(N)` was STALE.  The even-N guard existed because the old
//      m0_multiplier was wrong on multi-term inputs; M0MultiplierExact replaced it and
//      applies at EVERY level parity (SchoferFormula.m, the Nprimes branch).  The old line
//      mis-scored every ODD-N base, reporting a true multiplier of 0 + res instead of
//      applied + res.
//
//  (2) the verdict `mm eq tot` was VACUOUS.  With tot := mm + res and res = 0 it compares
//      mm to itself and always prints OK.  It also used M0Multiplier, the superseded
//      Route-C form that is correct on 15_2 only, rather than the M0MultiplierExact the
//      pipeline actually applies.
//
// The real test is res = 0: since res = true - applied and `applied` IS the value the
// pipeline used, res = 0 says exactly that the shipped multiplier is right on this base.
Ld := ShimuraCurveLattice(D, N);
applied := (N gt 1);
printf "the pipeline applies m0_multiplier here: %o  (every parity)\n", applied;
mexact := applied select M0MultiplierExact([fs[k] : k in [-1,-2]], Ld, D, N)
                 else [Rationals()!0, Rationals()!0];
for i -> k in [-1, -2] do
    ok := res[k] eq 0;
    printf "  form %-3o: residual %-6o applied %-8o TRUE multiplier %-8o %o\n",
           k, res[k], mexact[i], mexact[i] + res[k],
           ok select "pipeline OK" else "*** PIPELINE WRONG: true = applied + residual ***";
end for;
printf "TOTAL %o s\n", Cputime(t0);
quit;
