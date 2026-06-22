// Focused analysis of the OPEN (undetermined) cases in the completed run.
// Writes data/open_cases_analysis.txt
SetQuitOnError(true);
AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 0);

curves := eval Read("data/curves_after_UpdateCurves8.dat");
out := "data/open_cases_analysis.txt";
procedure W(s) Write(out, s); end procedure;

function Status(X)
    if not assigned X`IsSubhyp then return "open"; end if;
    return X`IsSubhyp select "proved" else "ruled";
end function;

// id -> curve index, for cover-lattice lookups
byid := AssociativeArray();
for i in [1..#curves] do byid[curves[i]`CurveID] := i; end for;

open_all := [X : X in curves | Status(X) eq "open"];
// Final step: D=1 quotients are classical modular curves X_0(N) whose
// (sub)hyperellipticity is already known in the literature (Ogg et al.),
// so they are resolved and excluded from the genuinely-open count.
open_D1 := [X : X in open_all | X`D eq 1];
open := [X : X in open_all | X`D ne 1];

W("================================================================");
W("OPEN (UNDETERMINED) CASES — ANALYSIS");
W("================================================================");
W(Sprintf("open before final step                 : %o", #open_all));
W(Sprintf("  resolved: D=1 modular curves known   : %o", #open_D1));
W(Sprintf("GENUINELY OPEN (D>1 quaternion Shimura) : %o", #open));
W("================================================================");
W("");

// ---- by genus ----
W("By genus");
W("--------");
for g in Sort(SetToSequence({X`g : X in open})) do
    W(Sprintf("  genus %o : %o", g, #[X:X in open|X`g eq g]));
end for;
W("");

// ---- by |W| (size of the AL group quotiented out) ----
W("By |W| (number of AL involutions in the quotient group, incl. identity)");
W("----------------------------------------------------------------------");
for w in Sort(SetToSequence({#X`W : X in open})) do
    W(Sprintf("  |W|=%-3o : %o", w, #[X:X in open|#X`W eq w]));
end for;
W("");

// ---- by number of AL involutions still acting non-trivially on the quotient ----
// The full AL group has order 2^omega(DN); quotienting by W leaves a group of order
// 2^omega(DN)/|W| acting on X/W, i.e. 2^omega(DN)/|W| - 1 non-trivial involutions.
W("By # AL involutions still acting non-trivially on the quotient (2^omega(DN)/|W| - 1)");
W("-----------------------------------------------------------------------------------");
actmult := AssociativeArray();
for X in open do
    k := (2^#PrimeDivisors(X`D * X`N) div #X`W) - 1;
    if not IsDefined(actmult, k) then actmult[k] := 0; end if;
    actmult[k] +:= 1;
end for;
actkeys := Sort([k : k in Keys(actmult)]);
W(Sprintf("  multiset: { %o }",
    Join([Sprintf("%o^%o", k, actmult[k]) : k in actkeys], ", ")));
for k in actkeys do
    W(Sprintf("  %o acting : %o", k, actmult[k]));
end for;
W("");

// ---- by D (indefinite quaternion discriminant) ----
W("By D (indefinite quaternion discriminant)");
W("-----------------------------------------");
for D in Sort(SetToSequence({X`D : X in open})) do
    W(Sprintf("  D=%-6o : %o", D, #[X:X in open|X`D eq D]));
end for;
W("");

// ---- distinct (D,N) families represented ----
fams := {<X`D, X`N> : X in open};
W(Sprintf("Distinct (D,N) families with at least one open quotient : %o", #fams));
famseq := Sort(SetToSequence(fams), func<a,b| a[1] ne b[1] select a[1]-b[1] else a[2]-b[2]>);
W("  (D, N, #open, total quotients of this (D,N), genus of X_0(D,N))");
for f in famseq do
    nopen := #[X:X in open|X`D eq f[1] and X`N eq f[2]];
    ntot  := #[X:X in curves|X`D eq f[1] and X`N eq f[2]];
    W(Sprintf("    D=%-5o N=%-5o  open=%-3o / %-3o quotients", f[1], f[2], nopen, ntot));
end for;
W("");

// ---- the genus-3 question: hyperelliptic vs smooth plane quartic ----
g3 := [X:X in open|X`g eq 3];
W("================================================================");
W(Sprintf("GENUS-3 open cases : %o", #g3));
W("A genus-3 curve is either hyperelliptic or a smooth plane quartic.");
W("These are open because no filter settled which; resolving them needs");
W("an explicit model / rational-point or gonality argument.");
W("================================================================");
W(Sprintf("  %-5o %-6o %-30o", "D", "N", "W"));
g3s := Sort(g3, func<a,b| a`D ne b`D select a`D-b`D else (a`N ne b`N select a`N-b`N else 0)>);
for X in g3s do
    W(Sprintf("  %-5o %-6o %-30o", X`D, X`N, Sort(SetToSequence(X`W))));
end for;
W("");

// ---- cover-lattice context: status of neighbours of each open curve ----
// Covers   = ids of curves THIS one maps ONTO (further quotients, smaller)
// CoveredBy = ids of curves that map ONTO this one (less-quotiented, larger)
W("================================================================");
W("COVER-LATTICE CONTEXT OF OPEN CURVES");
W("For each open curve, look at the curves it covers (its further");
W("quotients) and the curves that cover it (its parents).");
W("================================================================");
function NbStatus(ids)
    s := AssociativeArray();
    s["open"]:=0; s["proved"]:=0; s["ruled"]:=0;
    for id in ids do
        if IsDefined(byid, id) then
            st := Status(curves[byid[id]]);
            s[st] +:= 1;
        end if;
    end for;
    return s;
end function;

// How many open curves have a proved-subhyp curve among the things they cover?
covers_a_proved := 0; covered_by_a_proved := 0;
covers_only_ruled := 0; isolated := 0;
for X in open do
    cv  := assigned X`Covers   select X`Covers   else [];
    cb  := assigned X`CoveredBy select X`CoveredBy else [];
    sc := NbStatus(cv);
    sp := NbStatus(cb);
    if sc["proved"] gt 0 then covers_a_proved +:= 1; end if;
    if sp["proved"] gt 0 then covered_by_a_proved +:= 1; end if;
    if #cv gt 0 and sc["ruled"] eq #[i:i in cv|IsDefined(byid,i)] and sc["ruled"] gt 0 then
        covers_only_ruled +:= 1;
    end if;
    if #cv eq 0 and #cb eq 0 then isolated +:= 1; end if;
end for;
W(Sprintf("  open curves that COVER a proved sub-hyperelliptic quotient : %o", covers_a_proved));
W(Sprintf("  open curves COVERED BY a proved sub-hyperelliptic quotient : %o", covered_by_a_proved));
W(Sprintf("  open curves all of whose further quotients are ruled out   : %o", covers_only_ruled));
W(Sprintf("  open curves with no recorded cover relations (isolated)    : %o", isolated));
W("");

// ---- relations of open curves to OTHER open curves ----
// How many open curves could still be resolved "for free" by relating them to
// another open curve, rather than by an individual computation:
//   (a) it covers another open curve            -> upward closure may apply;
//   (b) it is isomorphic over F_p to an open curve (SpecialFiberIsomorphism: the
//       bigger level X_0(D,Np)/W whose normalization mod p is the open source
//       X_0(D,N)/W', W' = the AL in W prime to p);
//   (c) it is fully isomorphic (FH99 Prop 4/5) to an open curve.
// Counts are of distinct open curves, disjoint, taken in the order (a), (b), (c).
W("================================================================");
W("RELATIONS OF OPEN CURVES TO OTHER OPEN CURVES");
W("================================================================");
lut := AssociativeArray();
for i in [1..#curves] do c := curves[i]; lut[<c`D, c`N, c`W>] := i; end for;
openids := {X`CurveID : X in open};

coversOpen := {X`CurveID : X in open
                | assigned X`Covers and exists{id : id in X`Covers | id in openids}};

fpIso := {};
for X in open do
    for p in PrimeDivisors(X`N) do
        if Valuation(X`N, p) ne 1 then continue; end if;
        if not exists{w : w in X`W | w mod p eq 0} then continue; end if;
        key := <X`D, X`N div p, {w : w in X`W | w mod p ne 0}>;
        if IsDefined(lut, key) and curves[lut[key]]`CurveID in openids then
            Include(~fpIso, X`CurveID); break;
        end if;
    end for;
end for;

fullIso := {};
for X in open do
    for finder in [FindIsomorphicCurveProp4, FindIsomorphicCurveProp5] do
        Y := finder(X, curves, lut);
        if Type(Y) eq ShimuraQuot and Y`CurveID in openids then
            Include(~fullIso, X`CurveID);
        end if;
    end for;
end for;

bIso := fpIso diff coversOpen;
cIso := (fullIso diff coversOpen) diff fpIso;
W(Sprintf("  cover another open curve (closure may apply)         : %o", #coversOpen));
W(Sprintf("  + isomorphic over F_p to another open curve          : %o", #bIso));
W(Sprintf("  + fully isomorphic (FH99 Prop 4/5) to another open   : %o", #cIso));
W("");
W(Sprintf("  open curves with NO such lifeline -> must be resolved individually : %o",
    #open - #(coversOpen join fpIso join fullIso)));
W("");

printf "wrote %o\n", out;
quit;
