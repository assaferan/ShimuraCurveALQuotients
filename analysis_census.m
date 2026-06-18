// One-off census of the completed pipeline output.
// Writes a human-readable report to data/pipeline_analysis.txt.
SetQuitOnError(true);
AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 0);

curves := eval Read("data/curves_after_UpdateCurves8.dat");
out := "data/pipeline_analysis_census.txt";

procedure W(s)
    Write(out, s);
end procedure;

// ---- classification helpers ----
function Status(X)
    if not assigned X`IsSubhyp then return "open"; end if;
    if X`IsSubhyp then return "proved"; else return "ruled"; end if;
end function;

function Kind(X)
    if assigned X`IsP1 and X`IsP1 then return "P1"; end if;
    if assigned X`IsEC and X`IsEC then return "elliptic"; end if;
    if assigned X`IsHyp and X`IsHyp then return "hyperelliptic"; end if;
    return "?";
end function;

n := #curves;
open := [X : X in curves | Status(X) eq "open"];
ruled := [X : X in curves | Status(X) eq "ruled"];
proved := [X : X in curves | Status(X) eq "proved"];

W(Sprintf("FINAL CENSUS  (data/curves_after_UpdateCurves8.dat)"));
W(Sprintf("================================================================"));
W(Sprintf("Total AL quotients         : %o", n));
W(Sprintf("  proved sub-hyperelliptic : %o", #proved));
W(Sprintf("  ruled out (not subhyp)   : %o", #ruled));
W(Sprintf("  open (undetermined)      : %o", #open));
W("");

// ---- proved breakdown ----
p1   := [X : X in proved | Kind(X) eq "P1"];
ec   := [X : X in proved | Kind(X) eq "elliptic"];
hyp  := [X : X in proved | Kind(X) eq "hyperelliptic"];
W("Proved sub-hyperelliptic, by type");
W("---------------------------------");
W(Sprintf("  P^1 (genus 0)         : %o", #p1));
W(Sprintf("  elliptic (genus 1)    : %o", #ec));
W(Sprintf("  hyperelliptic (g>=2)  : %o", #hyp));
W("");

// ---- hyperelliptic by genus ----
W("Hyperelliptic quotients, by genus");
W("---------------------------------");
genera := Sort(SetToSequence({X`g : X in hyp}));
for g in genera do
    W(Sprintf("  genus %o : %o", g, #[X : X in hyp | X`g eq g]));
end for;
W("");

// ---- genus distribution over ALL curves ----
W(Sprintf("Genus distribution over ALL %o quotients (proved / ruled / open)", n));
allg := Sort(SetToSequence({X`g : X in curves}));
W(Sprintf("  %-6o %-8o %-8o %-8o %-8o", "genus", "total", "proved", "ruled", "open"));
for g in allg do
    cg := [X : X in curves | X`g eq g];
    W(Sprintf("  %-6o %-8o %-8o %-8o %-8o", g, #cg,
        #[X:X in cg|Status(X) eq "proved"],
        #[X:X in cg|Status(X) eq "ruled"],
        #[X:X in cg|Status(X) eq "open"]));
end for;
W("");

// ---- attribution: which test proved/ruled each curve ----
function TestTag(X)
    if not assigned X`TestInWhichProved then return "(none recorded)"; end if;
    t := X`TestInWhichProved;
    // group by the leading word (drops trailing curve-ids / involution names / counts)
    toks := Split(t, " ");
    return #toks eq 0 select "(none recorded)" else toks[1];
end function;

W("Attribution: which test recorded the determination (ruled-out curves)");
W("--------------------------------------------------------------------");
tags := {};
for X in ruled do Include(~tags, TestTag(X)); end for;
counts := [<t, #[X:X in ruled|TestTag(X) eq t]> : t in tags];
Sort(~counts, func<a,b | b[2]-a[2]>);
for c in counts do
    W(Sprintf("  %-45o : %o", c[1], c[2]));
end for;
W("");

W("Attribution: which test recorded the determination (proved subhyp)");
W("------------------------------------------------------------------");
tags := {};
for X in proved do Include(~tags, TestTag(X)); end for;
counts := [<t, #[X:X in proved|TestTag(X) eq t]> : t in tags];
Sort(~counts, func<a,b | b[2]-a[2]>);
for c in counts do
    W(Sprintf("  %-45o : %o", c[1], c[2]));
end for;
W("");

// ---- the open cases ----
W("OPEN (undetermined) cases — full list");
W("=====================================");
W(Sprintf("count = %o", #open));
W(Sprintf("  %-6o %-6o %-6o %-30o", "D", "N", "g", "W"));
opensort := open;
Sort(~opensort, func<a,b | a`D ne b`D select a`D-b`D else (a`N ne b`N select a`N-b`N else a`g-b`g)>);
for X in opensort do
    W(Sprintf("  %-6o %-6o %-6o %-30o", X`D, X`N, X`g, Sort(SetToSequence(X`W))));
end for;
W("");

W("Open cases, by genus");
W("--------------------");
for g in Sort(SetToSequence({X`g : X in open})) do
    W(Sprintf("  genus %o : %o", g, #[X:X in open|X`g eq g]));
end for;
W("");

W("Open cases, by D (quaternion discriminant)");
W("------------------------------------------");
for D in Sort(SetToSequence({X`D : X in open})) do
    W(Sprintf("  D=%-6o : %o", D, #[X:X in open|X`D eq D]));
end for;

printf "wrote %o\n", out;
quit;
