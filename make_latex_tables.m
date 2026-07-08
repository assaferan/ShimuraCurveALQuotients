// Auto-generate the two LaTeX filtering tables (star curves; all AL quotients)
// for the paper, directly from the pipeline output files data/curves_after_*.dat.
// Writes data/pipeline_tables.tex.  Run after the pipeline (or via
// make_pipeline_summary.sh).  Columns: Stage | Open | Filt. | Proved |
// New filt. | New proved, where the cumulative (Open,Filt,Proved) counts are read
// from the .dat snapshot at the end of each (group of) stage(s), and the "New"
// columns are the deltas from the previous row.
//
// Rows are grouped BY PAPER SECTION, not one-per-pipeline-stage:
//   * "Special fiber isomorphisms" merges HHProposition1 + SpecialFiberIsomorphism.
//   * "Weil polynomials" (star) merges FilterByWeilPolynomialStar + the redundant
//     FilterStarCurvesByFpAutomorphisms cross-check.
//   * "Refined Atkin--Lehner fixed points" merges FilterByComplicatedAL... +
//     FilterByGeneralizedComplicatedFixedPoints (both Section ref:fixedpointsAL).
// The three expansion operations (generate / genus<=2 / AL-as-hyperelliptic) are
// logical sub-steps of the single UpdateByGenus stage, so their intermediate
// counts are computed rather than read from a snapshot.  The final
// "Resolve classical modular curve cases from literature" row moves every still-open
// D=1 quotient (a classical modular curve, resolved by the literature) to Filt.
SetQuitOnError(true);
AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 0);

out := "data/pipeline_tables.tex";

// cumulative open / ruled(=Filt) / proved / total in a snapshot file
function CO(f)
    cs := eval Read(f);
    o := #[X : X in cs | not assigned X`IsSubhyp];
    p := #[X : X in cs | assigned X`IsSubhyp and X`IsSubhyp];
    r := #cs - o - p;
    return o, r, p, #cs;
end function;

// format one table row: label & open & filt & proved & newfilt & newproved
function Row(label, o, r, p, nf, np)
    return Sprintf("%o\n& %o & %o & %o & %o & %o \\\\", label, o, r, p, nf, np);
end function;

// ---------------------------------------------------------------------------
// TABLE 1 : star curves
// ---------------------------------------------------------------------------
T1 := [
 <"Reduction to finitely many pairs $(D,N)$ (Section~\\ref{sec:reduction})",
        "data/curves_after_UpdateGenera.dat">,
 <"Genus $\\leq 2$ curves (Section~\\ref{sec:genus})",
        "data/curves_after_UpdateByGenusStar.dat">,
 <"Finite field point count (Section~\\ref{sec:trace})",
        "data/curves_after_FilterByTraceStar.dat">,
 <"Special fiber isomorphisms (Section~\\ref{sec:specialfiber})",
        "data/curves_after_SpecialFiberIsomorphismStar.dat">,
 <"Weil polynomials (Section~\\ref{sec:Weilpolys})",
        "data/curves_after_FilterStarCurvesByFpAutomorphisms.dat">,
 <"Modular non-Atkin--Lehner involutions (Section~\\ref{sec:modularnonAL})",
        "data/curves_after_FilterByNonALInvolutionsStar.dat">
];

tab1 := "";
prevr := 0; prevp := 0; ntot := 0;
for i in [1..#T1] do
    o, r, p, n := CO(T1[i][2]);
    ntot := n;
    if i eq 1 then nf := "--"; np := "--";
    else nf := Sprint(r-prevr); np := Sprint(p-prevp); end if;
    tab1 cat:= Row(T1[i][1], o, r, p, nf, np) cat "\n";
    prevr := r; prevp := p;
end for;
nstar := ntot;

// ---------------------------------------------------------------------------
// TABLE 2 : all AL quotients (expansion folded in; no separate EXPANSION section)
// ---------------------------------------------------------------------------
// Phase-A end state (carried into the expanded set as inherited star verdicts).
provedA := 0; ruledA := 0; genus2stars := 0;
stars := eval Read("data/curves_after_FilterByNonALInvolutionsStar.dat");
provedA := #[X : X in stars | assigned X`IsSubhyp and X`IsSubhyp];
ruledA  := #stars - #[X : X in stars | not assigned X`IsSubhyp] - provedA;
genus2stars := #[X : X in stars | X`g le 2];

// Final expanded set (UpdateByGenus) = state after all three expansion ops.
o3, r3, p3, total := CO("data/curves_after_UpdateByGenus.dat");
expanded := eval Read("data/curves_after_UpdateByGenus.dat");
G2 := #[X : X in expanded | X`g le 2];            // all genus<=2 quotients
Pgenus := provedA + (G2 - genus2stars);            // proved after genus reclassification

tab2 := "";
// expansion sub-rows (computed)
tab2 cat:= Row("Generate all quotients $\\Curve{D}{N}/W$",
               total-ruledA-provedA, ruledA, provedA, "--", "--") cat "\n";
tab2 cat:= Row("Genus $\\leq 2$ curves (Section~\\ref{sec:genus})",
               total-ruledA-Pgenus, ruledA, Pgenus, 0, Pgenus-provedA) cat "\n";
tab2 cat:= Row("Atkin--Lehner involutions as hyperelliptic involutions (Section~\\ref{sec:ALashyp})",
               o3, r3, p3, r3-ruledA, p3-Pgenus) cat "\n";

// remaining Phase-B snapshots
T2 := [
 <"Propagate closure and isomorphism", "data/curves_after_UpdateCurves1.dat">,
 <"Atkin--Lehner fixed points (Section~\\ref{sec:fixedpointsAL})",
        "data/curves_after_FilterByALFixedPointsOnQuotient.dat">,
 <"Propagate closure and isomorphism", "data/curves_after_UpdateCurves2.dat">,
 <"Genus $3$ covers of genus $2$ curves (Section~\\ref{sec:g3coverg2})",
        "data/curves_after_Genus3CoversGenus2.dat">,
 <"Propagate closure and isomorphism", "data/curves_after_UpdateCurves3.dat">,
 <"Degeneracy morphism (Section~\\ref{sec:degeneracy})",
        "data/curves_after_FilterByDegeneracyMorphism.dat">,
 <"Propagate closure and isomorphism", "data/curves_after_UpdateCurves4.dat">,
 <"Refined Atkin--Lehner fixed points (Section~\\ref{sec:fixedpointsAL})",
        "data/curves_after_FilterByGeneralizedComplicatedFixedPoints.dat">,
 <"Propagate closure and isomorphism", "data/curves_after_UpdateCurves5.dat">,
 <"Finite field point count (Section~\\ref{sec:trace})",
        "data/curves_after_FilterByTrace.dat">,
 <"Propagate closure and isomorphism", "data/curves_after_UpdateCurves6.dat">,
 <"Weil polynomials (Section~\\ref{sec:Weilpolys})",
        "data/curves_after_FilterByWeilPolynomial.dat">,
 <"Propagate closure and isomorphism", "data/curves_after_UpdateCurves7.dat">,
 <"Modular non-Atkin--Lehner involutions (Section~\\ref{sec:modularnonAL})",
        "data/curves_after_FilterByNonALInvolutions.dat">,
 <"Propagate closure and isomorphism", "data/curves_after_UpdateCurves8.dat">
];

prevr := r3; prevp := p3; oF := 0; rF := 0; pF := 0;
for i in [1..#T2] do
    o, r, p, _ := CO(T2[i][2]);
    tab2 cat:= Row(T2[i][1], o, r, p, r-prevr, p-prevp) cat "\n";
    prevr := r; prevp := p;
    oF := o; rF := r; pF := p;
end for;

// literature row: every remaining D=1 quotient (classical modular curve) is
// resolved (non-sub-hyperelliptic) by the literature and moves Open -> Filt.
final := eval Read("data/curves_after_UpdateCurves8.dat");
d1open := #[X : X in final | not assigned X`IsSubhyp and X`D eq 1];
tab2 cat:= Row("Resolve classical modular curve cases from literature",
               oF-d1open, rF+d1open, pF, d1open, 0) cat "\n";
nquot := total;

// ---------------------------------------------------------------------------
// emit
// ---------------------------------------------------------------------------
pre1 := "\\begin{table}[ht]\n\\centering\n\\small\n \\renewcommand{\\arraystretch}{1.15}\n  \\setlength{\\tabcolsep}{3pt}\n\\begin{tabularx}{\\textwidth}{Xrrrrr}\n\\hline\nStage & Open & Filt. & Proved & New filt. & New proved \\\\\n\\hline";
post1 := Sprintf("\\hline\n\\end{tabularx}\n\\caption{Filtering the $%o$ star curves $\\StarCurve{D}{N}$.}\n\\label{tab:star-curve-filtering}\n\\end{table}", nstar);

pre2 := "\\begin{table}[ht]\n\\centering \\small \n \\renewcommand{\\arraystretch}{1.15}\n  \\setlength{\\tabcolsep}{3pt}\n\\begin{tabularx}{\\textwidth}{Xrrrrr}\n\\toprule\nOperation & Open & Filt. & Proved & New filt. & New proved \\\\\n\\midrule";
post2 := Sprintf("\\hline\n\\end{tabularx}\n\\caption{Filtering all $%o$ Atkin--Lehner quotients after expansion.}\n\\label{tab:quotient-filtering}\n\\end{table}", nquot);

tex := Sprintf("%% Auto-generated by make_latex_tables.m -- do not edit by hand.\n\n%o\n%o%o\n\n%o\n%o%o\n",
               pre1, tab1, post1, pre2, tab2, post2);

Write(out, tex : Overwrite := true);
printf "wrote %o\n", out;
quit;
