// Worker for parallel filtering of Shimura curve quotients.
//
// Run via:
//   magma input_dat:=... chunk:=N total_chunks:=M output_dat:=... [stage:=FilterByTrace] parallel_filter_worker.m
//
// Note: Magma receives command-line :=  values as raw strings (MonStgElt),
// so do NOT wrap path/stage values in Magma string quotes when calling.
// Integer args are converted with StringToInteger() below.
//
// The worker loads the full curve list, selects its cost-balanced subset of curves
// (see the assignment below), runs the named filter on them, then writes the (possibly
// modified) curves tagged with their original indices to output_dat.  parallel_merge.m
// reassembles all chunks back into the original order by index.
//
// Safe stages (per-curve computation, no global index access):
//   FilterByTrace, FilterByTraceStar,
//   FilterByALFixedPointsOnQuotient,
//   FilterByComplicatedALFixedPointsOnQuotient,
//   FilterByGeneralizedComplicatedFixedPoints, FilterByGeneralizedComplicatedFixedPointsStar,
//   FilterByWeilPolynomial, FilterByWeilPolynomialStar,
//   FilterByDegeneracyMorphism
//
// FilterStarCurvesByFpAutomorphisms is also safe (uses loop index, not CurveID)
//
// FilterByAutomorphismGroup is per curve, like the above.  FilterByTwistedTrace and
// FilterByTwistedWeilPolynomial (and their *Star versions) compute the modular symbols of level D*N once per level and
// share them between the curves at that level, so for them the unit of work is a LEVEL: all
// curves with the same (D,N) go to the same chunk (see the assignment below).
//
// The non-star stages that have a *Star version skip the star curves once that version has run
// (see star_of below): it already ran the same function on them.

// Convert integer args from the raw strings Magma receives on the command line
chunk_i        := StringToInteger(chunk);
total_chunks_i := StringToInteger(total_chunks);

if not assigned stage then
    stage := "FilterByTrace";
end if;

SetQuitOnError(true);
AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 0);

try

curves := eval Read(input_dat);
n := #curves;

// A non-star stage calls the same function as its star version, and a star curve's (D,N,W,g) do
// not change between the two, so re-running it on a star curve the star version left undecided
// repeats that computation and returns the same answer.  Skip those curves, but only when the star
// version has run in this data dir: its output curves_after_<stage>Star.dat sits next to
// input_dat, and the curve's level (D,N) is in it.
star_of := AssociativeArray();
star_of["FilterByTrace"]                 := "FilterByTraceStar";
star_of["FilterByTwistedTrace"]          := "FilterByTwistedTraceStar";
star_of["FilterByWeilPolynomial"]        := "FilterByWeilPolynomialStar";
star_of["FilterByTwistedWeilPolynomial"] := "FilterByTwistedWeilPolynomialStar";
star_of["FilterByNonALInvolutions"]      := "FilterByNonALInvolutionsStar";
star_of["FilterByGeneralizedComplicatedFixedPoints"] := "FilterByGeneralizedComplicatedFixedPointsStar";
skip := [false : i in [1..n]];
if IsDefined(star_of, stage) then
    parts := Split(input_dat, "/");
    dir := #parts gt 1 select &cat[p cat "/" : p in parts[1..#parts-1]] else "";
    if input_dat[1] eq "/" then dir := "/" cat dir; end if;
    star_dat := dir cat "curves_after_" cat star_of[stage] cat ".dat";
    ok, _ := OpenTest(star_dat, "r");
    if ok then
        star_levels := {<X`D, X`N> : X in eval Read(star_dat)};
        skip := [#X`W eq 2^#PrimeDivisors(X`D*X`N) and <X`D, X`N> in star_levels : X in curves];
    end if;
end if;

// Cost-aware assignment: order all curves by descending cost estimate (CurveCostProxy),
// then deal them round-robin into total_chunks groups.  This (a) spreads the heavy curves
// across distinct chunks so no chunk gets several of them, and (b) puts the heaviest curves
// in the lowest-numbered chunks, which GNU parallel dispatches first.  Each worker computes
// the same ordering deterministically, so the chunks partition the curves with no overlap.
proxy := [CurveCostProxy(curves[i], stage) : i in [1..n]];
for i in [1..n] do if skip[i] then proxy[i] := 0; end if; end for;   // skipped: no cost
if stage in {"FilterByTwistedTrace", "FilterByTwistedWeilPolynomial",
             "FilterByTwistedTraceStar", "FilterByTwistedWeilPolynomialStar"} then
    // Level-grouped stages: the same cost-aware strided deal, over levels instead of curves.  A
    // level's cost is its most expensive curve (the modular symbols dominate and are shared), and
    // every curve of the level, decided or not, goes with it, so the chunks still partition 1..n.
    lv := AssociativeArray();
    for i in [1..n] do
        key := <curves[i]`D, curves[i]`N>;
        if not IsDefined(lv, key) then lv[key] := []; end if;
        Append(~lv[key], i);
    end for;
    groups := [lv[key] : key in Keys(lv)];
    gcost := [Max([proxy[i] : i in grp]) : grp in groups];
    gfirst := [Min(grp) : grp in groups];
    gperm := [1..#groups];
    Sort(~gperm, func<a, b | gcost[a] gt gcost[b] select -1 else (gcost[a] lt gcost[b] select 1 else gfirst[a] - gfirst[b])>);
    my_idx := Sort(&cat([groups[gperm[k]] : k in [chunk_i .. #groups by total_chunks_i]] cat [[Integers()|]]));
else
    perm := [1..n];
    Sort(~perm, func<i, j | proxy[i] gt proxy[j] select -1 else (proxy[i] lt proxy[j] select 1 else i - j)>);
    my_idx := [perm[k] : k in [chunk_i .. n by total_chunks_i]];   // strided slice of the sorted order
end if;
subseq := [curves[i] : i in my_idx];
// Hand the filter only the curves it is not skipping; the skipped ones are written back unchanged.
full_subseq := subseq;
keep := [j : j in [1..#subseq] | not skip[my_idx[j]]];
subseq := [full_subseq[j] : j in keep];

t0 := Realtime();

case stage:
    when "FilterByTrace", "FilterByTraceStar":
        FilterByTrace(~subseq);
    when "FilterByAutomorphismGroup":
        FilterByAutomorphismGroup(~subseq);
    when "FilterByTwistedTrace", "FilterByTwistedTraceStar":
        FilterByTwistedTrace(~subseq);
    when "FilterByTwistedWeilPolynomial", "FilterByTwistedWeilPolynomialStar":
        FilterByTwistedWeilPolynomial(~subseq);
    when "FilterStarCurvesByFpAutomorphisms":
        FilterStarCurvesByFpAutomorphisms(~subseq);
    when "FilterByALFixedPointsOnQuotient":
        FilterByALFixedPointsOnQuotient(~subseq);
    when "FilterByComplicatedALFixedPointsOnQuotient":
        FilterByComplicatedALFixedPointsOnQuotient(~subseq);
    when "FilterByGeneralizedComplicatedFixedPoints", "FilterByGeneralizedComplicatedFixedPointsStar":
        FilterByGeneralizedComplicatedFixedPoints(~subseq);
    when "FilterBySpecialFiber":
        FilterBySpecialFiber(~subseq);
    when "FilterByDegeneracyMorphism":
        FilterByDegeneracyMorphism(~subseq);
    when "FilterByWeilPolynomial", "FilterByWeilPolynomialStar":
        FilterByWeilPolynomialGenusScaled(~subseq);
    when "FilterByNonALInvolutions", "FilterByNonALInvolutionsStar":
        FilterByNonALInvolutions(~subseq);
    else
        error Sprintf("Unknown or unsupported stage: %o", stage);
end case;

catch e
    WriteStderr(Sprintf("ERROR in worker stage %o chunk %o/%o:\n", stage, chunk_i, total_chunks_i));
    WriteStderr(e);
    error e;  // re-raise so SetQuitOnError exits non-zero
end try;

for t->j in keep do full_subseq[j] := subseq[t]; end for;
subseq := full_subseq;

// Tag each curve with its original index so the (index-aware) merge can restore order.
Write(output_dat, Sprint([<my_idx[j], subseq[j]> : j in [1..#subseq]], "Magma") : Overwrite);
printf "Worker %o/%o: %o curves (%o star curves skipped), %o s\n", chunk_i, total_chunks_i, #subseq, #subseq - #keep, Realtime() - t0;
quit;
