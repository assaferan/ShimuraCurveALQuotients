// Merges the index-tagged chunk files produced by parallel_filter_worker.m back into a
// single .dat file in the original curve order.
//
// Usage:
//   magma chunks_dir:=... total_chunks:=... output_dat:=... parallel_merge.m
//
// Each chunk file is a Magma sequence of <original_index, ShimuraQuot> tuples:
//   [ <i1, CreateShimuraQuot(...)>, <i2, CreateShimuraQuot(...)>, ... ]
// (since the workers now distribute curves cost-aware/strided rather than contiguously,
//  the merge must reorder by original index rather than just concatenate).

SetQuitOnError(true);
AttachSpec("ShimuraQuotients.spec");

tc := StringToInteger(total_chunks);

tagged := [];
for c in [1..tc] do
    path := chunks_dir cat "/chunk_" cat IntegerToString(c) cat "_of_" cat IntegerToString(tc) cat ".dat";
    part := eval Read(path);   // SeqEnum of <index, ShimuraQuot>
    tagged cat:= part;
end for;

Sort(~tagged, func<a, b | a[1] - b[1]>);   // restore original order
curves := [t[2] : t in tagged];

Write(output_dat, Sprint(curves, "Magma") : Overwrite);
printf "Merged %o curves into %o\n", #curves, output_dat;
quit;
