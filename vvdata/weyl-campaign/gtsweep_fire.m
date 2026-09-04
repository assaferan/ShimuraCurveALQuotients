// gtsweep_fire.m -- the FIRE variant of vvdata/gtsweep.m, kept HERE and not there.
//
// WHY IT LIVES ON THE CAMPAIGN BRANCH, IN THIS DIRECTORY. `vvdata/gtsweep.m` is a SHARED
// path: it exists on `main` too (tests/M0MultiplierExact.m cites it for ground-truth
// provenance, and main carries the gtsweep_*.log outputs). Shared-path files that diverge
// between main and this branch are exactly how two capability gaps went unnoticed --
// nmzsolve.py's t-shift fallback was missing from main for nine days. So the shared copy is
// kept byte-identical to main's, and this experiment lives in vvdata/weyl-campaign/, which
// main does not have and which therefore cannot diverge.
//
// WHAT IT IS. gtsweep.m plus the FIRE knob: admit up to NFIRE firing degree-1 discriminants
// into must_use, on the theory that NOT DETERMINED verdicts are starved of exactly those.
// Run it the same way, e.g.  magma -b DD:=34 NN:=3 FIRE:=3 vvdata/weyl-campaign/gtsweep_fire.m
//
// ⚠ THE THEORY WAS MEASURED AND IT IS WRONG -- see the block lower down. Kept because the
// knob is still the right PLACE to fix this, and because a measured negative is worth more
// than the confident comment it replaced. Do NOT port it to main on the old comment's word.
//
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

// FIRING DEGREE-1 DISCRIMINANTS -- what NOT DETERMINED actually lacks.
//
// A row is INFORMATIVE only when its firing status differs from a normaliser's, and
// HauptmodulM0Residuals discards every row with degs ne 1.  The Keep-admitted `extra` above
// are all NON-firing, so they add rows with exponents [0,0] that constrain nothing; and
// raising MaxNum does not help either, because the budget gets spent on degree-2 points
// (measured on X0^34(3): MN:=14 added six points, all degree 2, verdict unchanged).
//
// rat_all is the FIRST return of RationalandQuadraticCMPoints, i.e. the RATIONAL points, so
// everything in it is already degree 1.  Forcing a few FIRING ones into Include is therefore
// exactly the missing lever.  Three bases are stuck for want of them: X0^33(2) (one
// informative row, -15; its only other firing disc -55 is degree 2), X0^34(3) (two rows,
// both [1,1]) and X0^46(3) (r2 appears in no row at all).
//
// FIRE := 0 reproduces the earlier logs exactly.
//
// ⚠ MEASURED 2026-09-04, AND THE LEVER DOES NOT WORK. The claim above -- that these three
// bases are "stuck for want of" firing degree-1 discriminants and that admitting them is
// "exactly the missing lever" -- is NOT borne out. Ran each base at FIRE:=0 and FIRE:=3:
//
//     33_2   admitted [-15]        NOT MEASURABLE -> NOT MEASURABLE   (unchanged)
//     34_3   admitted [-11, -20]   NOT DETERMINED -> NOT DETERMINED   (unchanged)
//     46_3   admitted []           NOT DETERMINED, 25 pairs both ways (VACUOUS: its only
//                                  firing degree-1 candidate, -8, is already in must_use)
//
// So on the two bases where the knob actually admits anything, the verdict does not move;
// on the third it cannot engage at all -- consistent with the aside above that 46_3's r2
// "appears in no row at all", which is a different problem. The knob is kept (it is still
// the right PLACE to fix this) but it is UNVALIDATED -- do not cite it as a working lever,
// and do not port it to main on the strength of the comment above. This is also why main
// does not have it, while main DID get the t-shift fallback, which had real verification.
firing1 := [ p[1] : p in rat_all
             | not IsEmpty(PrimeDivisors(N div GCD(N, FundamentalDiscriminant(p[1])))) ];
NFIRE := 3; if assigned FIRE then NFIRE := StringToInteger(FIRE); end if;
extra_f := [ d : d in firing1 | d notin must_use ];
extra_f := extra_f[1 .. Minimum(NFIRE, #extra_f)];
printf "  firing degree-1 candidates %o; admitting %o\n", firing1, extra_f;
must_use := must_use join Set(extra_f);

// MN raises the point budget.  When the sweep returns NOT DETERMINED the missing
// information is FIRING discriminants: the Keep-admitted extras are all non-firing, so they
// contribute rows with exponents [0,0] that constrain nothing about (r1,r2), while the
// firing ones give [k-kA, k-kB] rows that do.  Firing discriminants arrive through the
// ordinary coprime-to-level channel and are capped by MaxNum, so raising it is the lever.
//
// MEASURED, AND IT DID NOT WORK: on X0^34(3), MN:=14 added six points and ALL SIX WERE
// DEGREE 2, which HauptmodulM0Residuals skips (`degs[i] ne 1`).  Informative rows stayed at
// two and the verdict stayed NOT DETERMINED.  Raising MaxNum alone is therefore not enough
// -- the budget has to be spent on DEGREE-1 firing discriminants, which this knob cannot
// ask for.  Kept because the knob is still the right place to fix that.
//
// It also RE-BASES the answer: the surviving pairs moved {<0,0>,<1,2>} -> {<-3,-3>,<-2,-1>},
// a uniform (-3,-3) shift, because the multiplier the pipeline APPLIED changed with the
// point set (non-firing table entries moved by 3^6, firing by 3^3).  Since the residual is
// true - applied, that is a change of origin and not a contradiction -- but it does mean a
// residual, and the "pipeline OK" verdict built on it, is relative to that run's point set.
MAXNUM := Maximum(7, #must_use + 2);
if assigned MN then MAXNUM := Maximum(MAXNUM, StringToInteger(MN)); end if;
all_cm_pts := CandidateDiscriminants(star, curves : Keep := must_use);
abs_tab, all_cm_pts := AbsoluteValuesAtCMPoints(star, curves, all_cm_pts, fs
                                                : MaxNum := MAXNUM,
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
