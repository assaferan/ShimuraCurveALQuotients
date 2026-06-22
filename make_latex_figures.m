// Auto-generate the two LaTeX figures for the paper from the pipeline output
// (data/curves_after_UpdateCurves8.dat).  Writes data/pipeline_figures.tex:
//   1. fig:genus-distribution  -- stacked bar chart of genus x status.
//   2. fig:test-attribution    -- per-test count of the resolved quotients.
//
// Status for the genus chart is the FINAL status: the still-open D=1 quotients
// are classical modular curves resolved (not sub-hyperelliptic) by the literature,
// so they are counted under "Proved not subhyperelliptic", leaving "Open" = the
// genuinely-open D>1 quotients.  The attribution figure lists, per resolving reason,
// every quotient that is not genuinely open: the pipeline TESTS (via TestInWhichProved)
// plus a "Classical modular curve (literature)" row for the still-open D=1 cases.
//
// Tests are grouped BY PAPER SECTION (same rule as make_latex_tables.m): the
// consecutive same-section stages are merged -- ComplicatedAL + Generalized ->
// "Refined Atkin--Lehner fixed points"; SpecialFiberIsomorphism + HHproposition1
// -> "Special fiber isomorphism".  Rows are sorted by count within each block and
// the log-scaled bars are normalised to the largest count.
SetQuitOnError(true);
AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 0);

curves := eval Read("data/curves_after_UpdateCurves8.dat");
out := "data/pipeline_figures.tex";
n := #curves;

function Status(X)
    if not assigned X`IsSubhyp then return "open"; end if;
    return X`IsSubhyp select "proved" else "ruled";
end function;

// ---------------------------------------------------------------------------
// FIGURE 1 : genus distribution (final status)
// ---------------------------------------------------------------------------
genera := Sort(SetToSequence({X`g : X in curves}));
subS := []; notS := []; opS := [];
for g in genera do
    cg := [X : X in curves | X`g eq g];
    sub := #[X : X in cg | Status(X) eq "proved"];
    rul := #[X : X in cg | Status(X) eq "ruled"];
    d1op := #[X : X in cg | Status(X) eq "open" and X`D eq 1];   // literature -> not subhyp
    op  := #[X : X in cg | Status(X) eq "open"] - d1op;
    Append(~subS, sub); Append(~notS, rul + d1op); Append(~opS, op);
end for;

xs := Join([Sprint(g) : g in genera], ",");
function Coords(seq)
    return &cat[Sprintf("(%o,%o) ", genera[i], seq[i]) : i in [1..#genera]];
end function;

fig1 :=
  "\\begin{figure}[ht]\n  \\centering\n  \\begin{tikzpicture}\n  \\begin{axis}[\n"
  cat "  ybar stacked,\n  width=\\textwidth,\n  height=0.50\\textwidth,\n  ymin=0,\n"
  cat "  xlabel={Genus},\n  ylabel={Curves},\n  symbolic x coords={" cat xs cat "},\n"
  cat "  xtick=data,\n  x tick label style={font=\\scriptsize},\n"
  cat "  legend style={at={(0.5,-0.20)}, anchor=north, legend columns=3},\n"
  cat "  title={Genus distribution of all $" cat Sprint(n) cat "$ quotients},\n  grid=major,\n  ]\n\n"
  cat "  \\addplot+[\n    draw=black,\n    line width=0.25pt,\n    fill=StatusSub!85\n  ] coordinates {\n  "
        cat Coords(subS) cat "\n  };\n\n"
  cat "  \\addplot+[\n    draw=black,\n    line width=0.25pt,\n    fill=StatusNot!65\n  ] coordinates {\n  "
        cat Coords(notS) cat "\n  };\n\n"
  cat "  \\addplot+[\n    draw=black,\n    line width=0.25pt,\n    fill=StatusOpen!45\n  ] coordinates {\n  "
        cat Coords(opS) cat "\n  };\n\n"
  cat "  \\legend{Proved subhyperelliptic, Proved not subhyperelliptic, Open}\n"
  cat "  \\end{axis}\n  \\end{tikzpicture}\n"
  cat "  \\caption{Final genus distribution of all quotients, separated by status.}\n"
  cat "  \\label{fig:genus-distribution}\n\\end{figure}\n";

// ---------------------------------------------------------------------------
// FIGURE 2 : test attribution
// ---------------------------------------------------------------------------
function Canon(X)
    if not assigned X`TestInWhichProved then return ""; end if;
    toks := Split(X`TestInWhichProved, " ");
    if #toks eq 0 then return ""; end if;
    t := toks[1];
    while #t gt 0 and (t[#t] eq "," or t[#t] eq ":") do t := t[1..#t-1]; end while;
    return t;
end function;

proved := [X : X in curves | Status(X) eq "proved"];
ruled  := [X : X in curves | Status(X) eq "ruled"];

provedGroups := [
 <"Hyperelliptic Atkin--Lehner involution", ["HyperellipticALInvolution"]>,
 <"Genus $\\leq 2$", ["GenusLt3"]>,
 <"Genus 3 cover of genus 2", ["Genus3CoverGenus2"]>,
 <"Modular non-Atkin--Lehner involution", ["ModularNonALInvolution"]>,
 <"Isomorphism", ["UpdateIsoStatus"]>
];
ruledGroups := [
 <"Atkin--Lehner fixed points", ["ALFixedPointsOnQuotient"]>,
 <"Finite field point count", ["Trace"]>,
 <"Refined Atkin--Lehner fixed points", ["ComplicatedALFixedPointsOnQuotient", "GeneralizedComplicatedFixedPoints"]>,
 <"Weil polynomial", ["WeilPolynomial"]>,
 <"Upward closure", ["UpwardClosure"]>,
 <"Special fiber isomorphism", ["SpecialFiberIsomorphism", "HHproposition1"]>,
 <"Isomorphism", ["UpdateIsoStatus"]>,
 <"Modular non-Atkin--Lehner involution", ["ModularNonALInvolution"]>,
 <"Degeneracy morphism", ["DegeneracyMorphism"]>
];

function Rows(curveset, groups)
    rows := [];
    for grp in groups do
        c := #[X : X in curveset | Canon(X) in grp[2]];
        if c gt 0 then Append(~rows, <grp[1], c>); end if;
    end for;
    Sort(~rows, func<a,b | b[2]-a[2]>);
    return rows;
end function;

provedRows := Rows(proved, provedGroups);
ruledRows  := Rows(ruled, ruledGroups);

// still-open D=1 quotients are classical modular curves resolved (not sub-
// hyperelliptic) by the literature -- record them as their own reason.
d1open := #[X : X in curves | Status(X) eq "open" and X`D eq 1];
if d1open gt 0 then
    Append(~ruledRows, <"Classical modular curve (literature)", d1open>);
    Sort(~ruledRows, func<a,b | b[2]-a[2]>);
end if;

maxcount := Max([r[2] : r in provedRows cat ruledRows]);
resolved := &+[r[2] : r in provedRows] + &+[r[2] : r in ruledRows];

macro :=
  "% \\testbar normalises log-scaled bars to the largest count; drop this line if\n"
  cat "% the macro is already defined in your preamble.\n"
  cat "\\newcommand{\\testbar}[2]{%\n"
  cat "  \\pgfmathsetmacro{\\barwidth}{max(0.6, 100*ln(#2)/ln(" cat Sprint(maxcount) cat "))}%\n"
  cat "  \\begin{tikzpicture}[baseline=-0.55ex]\n"
  cat "    \\fill[#1!10] (0,0) rectangle (2.7,0.16);\n"
  cat "    \\draw[black!45, line width=0.2pt] (0,0) rectangle (2.7,0.16);\n"
  cat "    \\fill[#1!75] (0,0) rectangle ({2.7*\\barwidth/100},0.16);\n"
  cat "  \\end{tikzpicture}%\n}\n";

function ARow(label, cnt, color)
    return label cat "\n  & " cat Sprint(cnt) cat " & \\testbar{" cat color cat "}{" cat Sprint(cnt) cat "} \\\\\n\n";
end function;

fig2 :=
  "\\begin{figure}[ht]\n\\centering\n\\small\n\\renewcommand{\\arraystretch}{1.12}\n\n"
  cat "\\begin{tabularx}{\\linewidth}{%\n  >{\\raggedright\\arraybackslash}p{0.47\\linewidth}\n  r\n"
  cat "  >{\\raggedright\\arraybackslash}p{0.26\\linewidth}\n}\n"
  cat "\\toprule\n\\textbf{Test} & \\textbf{Count} & \\textbf{Relative size} \\\\\n\\midrule\n"
  cat "\\multicolumn{3}{l}{\\textbf{Proved subhyperelliptic}} \\\\\n\\addlinespace[0.2em]\n\n";
for r in provedRows do fig2 cat:= ARow(r[1], r[2], "StatusSub"); end for;
fig2 cat:=
  "\\addlinespace[0.6em]\n\\midrule\n\\multicolumn{3}{l}{\\textbf{Proved not subhyperelliptic}} \\\\\n\\addlinespace[0.2em]\n\n";
for r in ruledRows do fig2 cat:= ARow(r[1], r[2], "StatusNot"); end for;
fig2 cat:=
  "\\bottomrule\n\\end{tabularx}\n\n"
  cat "\\caption{Distribution of the tests that determined the $" cat Sprint(resolved)
  cat "$ resolved quotients. Bar lengths are log-scaled.}\n\\label{fig:test-attribution}\n\\end{figure}\n";

// ---------------------------------------------------------------------------
tex := "% Auto-generated by make_latex_figures.m -- do not edit by hand.\n\n"
       cat fig1 cat "\n\n" cat macro cat "\n" cat fig2;
Write(out, tex : Overwrite := true);
printf "wrote %o  (proved %o, ruled %o, resolved-by-test %o, bar-max %o)\n",
       out, #proved, #ruled, resolved, maxcount;
quit;
