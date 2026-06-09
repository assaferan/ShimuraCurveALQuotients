// Worker for parallel filtering of Shimura curve quotients.
//
// Run via:
//   magma input_dat:=... chunk:=N total_chunks:=M output_dat:=... [stage:=FilterByTrace] parallel_filter_worker.m
//
// Note: Magma receives command-line :=  values as raw strings (MonStgElt),
// so do NOT wrap path/stage values in Magma string quotes when calling.
// Integer args are converted with StringToInteger() below.
//
// The worker loads the full curve list, extracts its assigned contiguous slice
// curves[start..end], runs the named filter on that slice, then writes the
// (possibly modified) slice to output_dat.  The merge script reassembles the
// slices in order.
//
// Safe stages (per-curve computation, no global index access):
//   FilterByTrace, FilterByTraceStar,
//   FilterByALFixedPointsOnQuotient,
//   FilterByComplicatedALFixedPointsOnQuotient,
//   FilterByWeilPolynomial,
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

chunk_size := (n + total_chunks_i - 1) div total_chunks_i;
start_idx  := (chunk_i - 1) * chunk_size + 1;
end_idx    := Minimum(chunk_i * chunk_size, n);

subseq := curves[start_idx .. end_idx];

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
    when "FilterByWeilPolynomial":
        FilterByWeilPolynomialGenusScaled(~subseq);
    when "FilterByNonALInvolutions":
        FilterByNonALInvolutions(~subseq);
    else
        error Sprintf("Unknown or unsupported stage: %o", stage);
end case;

catch e
    WriteStderr(Sprintf("ERROR in worker stage %o chunk %o/%o:\n", stage, chunk_i, total_chunks_i));
    WriteStderr(e);
    error e;  // re-raise so SetQuitOnError exits non-zero
end try;

Write(output_dat, Sprint(subseq, "Magma") : Overwrite);
printf "Worker %o/%o [%o..%o]: %o curves, %o s\n", chunk_i, total_chunks_i, start_idx, end_idx, #subseq, Realtime() - t0;
quit;
