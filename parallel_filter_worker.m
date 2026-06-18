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
//   FilterByWeilPolynomial, FilterByWeilPolynomialStar,
//   FilterByDegeneracyMorphism
//
// FilterStarCurvesByFpAutomorphisms is also safe (uses loop index, not CurveID)

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

// Cost-aware assignment: order all curves by descending cost estimate (CurveCostProxy),
// then deal them round-robin into total_chunks groups.  This (a) spreads the heavy curves
// across distinct chunks so no chunk gets several of them, and (b) puts the heaviest curves
// in the lowest-numbered chunks, which GNU parallel dispatches first.  Each worker computes
// the same ordering deterministically, so the chunks partition the curves with no overlap.
proxy := [CurveCostProxy(curves[i], stage) : i in [1..n]];
perm := [1..n];
Sort(~perm, func<i, j | proxy[i] gt proxy[j] select -1 else (proxy[i] lt proxy[j] select 1 else i - j)>);
my_idx := [perm[k] : k in [chunk_i .. n by total_chunks_i]];   // strided slice of the sorted order
subseq := [curves[i] : i in my_idx];

t0 := Realtime();

case stage:
    when "FilterByTrace", "FilterByTraceStar":
        FilterByTrace(~subseq);
    when "FilterStarCurvesByFpAutomorphisms":
        FilterStarCurvesByFpAutomorphisms(~subseq);
    when "FilterByALFixedPointsOnQuotient":
        FilterByALFixedPointsOnQuotient(~subseq);
    when "FilterByComplicatedALFixedPointsOnQuotient":
        FilterByComplicatedALFixedPointsOnQuotient(~subseq);
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

// Tag each curve with its original index so the (index-aware) merge can restore order.
Write(output_dat, Sprint([<my_idx[j], subseq[j]> : j in [1..#subseq]], "Magma") : Overwrite);
printf "Worker %o/%o: %o curves, %o s\n", chunk_i, total_chunks_i, #subseq, Realtime() - t0;
quit;
