// The curves X_0^D(N)/W with D > 1 whose stored model in data/models/ is PROVED, so that
// FilterByCanonicalQuadrics may decide them from it.  Agreement of a model's genus and point counts
// with the trace formula does not prove it, and neither does a model whose CM values rest on the
// m = 0 term of Schofer's formula at N > 1; such curves are not listed and stay undecided until a
// certificate decides them.
//
// Entries are <D, N, W, argument>: W the full Atkin-Lehner group as a sequence of indices m, and
// argument a short reference to what proves the equation (a paper or a note under paper/, and the
// test that reruns it).  Read with eval (Read(file) cat "\nreturn proved_models;").
proved_models := [];
