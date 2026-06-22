// Reconstruct the correct TestInWhichProved attribution for the final pipeline
// state.  Several filters used to overwrite TestInWhichProved when they re-touched
// an already-determined curve (without ever flipping the IsSubhyp value), so the
// final data/curves_after_UpdateCurves8.dat -- and the figure built from it --
// credited the wrong test.  Here we replay the saved per-stage files in pipeline
// order and, for each curve, keep the attribution of the FIRST stage that assigned
// IsSubhyp.  IsSubhyp/IsHyp are unchanged (determinations never flip); only
// TestInWhichProved is corrected.  Result is written to
// data/curves_after_UpdateCurves8_fixed.dat.
SetQuitOnError(true);
SetColumns(0);  // no line wrapping, so the map file has one entry per line
AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 0);

function Load(name)
    return eval Read(Sprintf("data/curves_after_%o.dat", name));
end function;

// Star stages run first (on X_0^*(D,N) only); matched to the full list by (D,N,W).
star_stages := [
 "UpdateByGenusStar",
 "FilterByTraceStar",
 "HHProposition1",
 "SpecialFiberIsomorphismStar",
 "FilterStarCurvesByFpAutomorphisms",
 "FilterByNonALInvolutionsStar"
];
// Full stages run after expansion; matched by CurveID (== index).
full_stages := [
 "UpdateByGenus",
 "UpdateCurves1",
 "FilterByALFixedPointsOnQuotient",
 "UpdateCurves2",
 "Genus3CoversGenus2",
 "UpdateCurves3",
 "FilterByDegeneracyMorphism",
 "UpdateCurves4",
 "FilterByComplicatedALFixedPointsOnQuotient",
 "UpdateCurves5",
 "FilterByTrace",
 "UpdateCurves6",
 "FilterByWeilPolynomial",
 "UpdateCurves7",
 "FilterByNonALInvolutions",
 "UpdateCurves8"
];

final := Load("UpdateCurves8");
n := #final;
printf "final list: %o curves\n", n;

// attribution[i] = corrected TestInWhichProved string for final[i];
// the BoolElt false is a sentinel meaning "not yet assigned".
attrib := [* false : i in [1..n] *];
firstStage := [* "" : i in [1..n] *];

// (D,N,W) -> final index, for matching star curves.
key2idx := AssociativeArray();
for i->X in final do
    key2idx[<X`D, X`N, X`W>] := i;
end for;

procedure consider(~attrib, ~firstStage, idx, X, stagename)
    if Type(attrib[idx]) ne BoolElt then return; end if;  // already recorded
    if not assigned X`IsSubhyp then return; end if;
    if not assigned X`TestInWhichProved then return; end if;
    attrib[idx] := X`TestInWhichProved;
    firstStage[idx] := stagename;
end procedure;

for s in star_stages do
    cs := Load(s);
    for X in cs do
        key := <X`D, X`N, X`W>;
        if not IsDefined(key2idx, key) then continue; end if;
        idx := key2idx[key];
        consider(~attrib, ~firstStage, idx, X, s);
    end for;
end for;

for s in full_stages do
    cs := Load(s);
    assert #cs eq n;
    for i in [1..n] do
        consider(~attrib, ~firstStage, i, cs[i], s);
    end for;
end for;

// Report: changes vs. the (buggy) final attribution.
nchanged := 0; nassigned := 0;
changeBy := AssociativeArray();   // "oldCanon -> newCanon" : count
function Canon(t)
    toks := Split(t, " ");
    if #toks eq 0 then return ""; end if;
    u := toks[1];
    while #u gt 0 and (u[#u] eq "," or u[#u] eq ":") do u := u[1..#u-1]; end while;
    return u;
end function;

// Emit map "CurveID<TAB>corrected attribution" for every assigned curve.
mapstr := "";
for i in [1..n] do
    if not assigned final[i]`IsSubhyp then continue; end if;
    nassigned +:= 1;
    error if Type(attrib[i]) eq BoolElt,
        Sprintf("curve %o is assigned in final but no stage file recorded its attribution", i);
    oldt := assigned final[i]`TestInWhichProved select final[i]`TestInWhichProved else "";
    newt := attrib[i];
    mapstr cat:= Sprintf("%o\t%o\n", final[i]`CurveID, newt);
    if Canon(oldt) ne Canon(newt) then
        nchanged +:= 1;
        k := Canon(oldt) cat "  ==>  " cat Canon(newt);
        if IsDefined(changeBy, k) then changeBy[k] +:= 1; else changeBy[k] := 1; end if;
    end if;
end for;

printf "assigned (resolved) curves: %o\n", nassigned;
printf "curves whose attribution token changed: %o\n", nchanged;
print "--- changes (oldToken ==> newToken : count) ---";
ks := Sort([k : k in Keys(changeBy)], func<a,b | changeBy[b]-changeBy[a]>);
for k in ks do printf "  %-70o %o\n", k, changeBy[k]; end for;

Write("data/attribution_map.txt", mapstr : Overwrite);
printf "wrote data/attribution_map.txt (%o entries)\n", nassigned;
quit;
