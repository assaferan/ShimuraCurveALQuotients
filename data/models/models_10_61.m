// Subhyperelliptic cover models for X_0(10,61)*
//
// PROVENANCE. Produced 2026-09-06 on lovelace, DEFAULT flags (M0PROGRESS=1 is diagnostics only):
//     NORMALIZ_BIN=... magma -b D_s:=10 N_s:=61 OUTDIR:=... genmodels.m < /dev/null
// AllEquationsAboveCovers took 147629 s (41 h); 7 cover-keys, ALL populated.
//
// WHY THIS ONE MATTERS. 10_61 is one of the OBSTRUCTED bases -- the base whose gate-3 failure was
// root-caused as two poisoned cosets (wi = 2 and wi = 1221), 924 failures with zero near-misses,
// where k0 overflowed to +inf at one coset and collapsed to ~1e-10 at the other. That was a real
// upstream defect, not a tolerance; with it cleared the base runs to completion. This is the first
// model from that class.
//
// ⚠ NOT a Guo-Yang base -- X_0^10(61) has no published equation, so there is NO external oracle
// here. What validates it is ModelChecks alone: genus self-consistency, the Shimura-curve genus
// formula, Weil-polynomial divisibility, and Eichler-Selberg trace-formula point counts, the last
// of which is independent of the Borcherds/Schofer path that produced the model. That is a weaker
// bar than the published-equation bases clear, and it should be quoted as such.
//
// ⚠ CODE VERSION, AND A PROCESS HAZARD WORTH KNOWING. The run started 2026-09-04 from lovelace's
// clone, pre-dating both the vx fix and the q-expansion bootstrap. It has NOT been regenerated on
// current main; re-deriving costs ~41 h.
// ⚠⚠ That clone was `git pull`-ed to main ON 2026-09-06 WHILE THIS JOB (and seven others) WERE
// STILL RUNNING FROM IT. `AttachSpec` loads packages ON DEMAND, so a long job can compile a
// package whose source changed under it and end up mixing code versions. No harm is evident here
// -- the run was already deep in AllEquationsAboveCovers, and ModelChecks passes 63/63 -- but the
// output cannot be pinned to a single commit, and that is a defect in how it was produced, not a
// property of the base.
// ⇒ RULE: launch long runs from a COPIED tree (as ~/shimura/vxfix does), never from the clone you
// intend to update. Do not update a clone that has jobs running from it.
//
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,5,122,610]] := [* <1, P![ 1/320, -33/1600, 57/1280, -7/200, 1/100 ], P![]> *];
models[[Integers()|1,2,305,610]] := [* <2, P![ 0, -1/400, 21/1280, -237/6400, 993/25600, -5/256, 1/256 ], P![]> *];
models[[Integers()|1,10,61,610]] := [* <4, P![ 0, -1/32000, 1053/2560000, -579/256000, 34879/5120000, -63679/5120000, 590237/40960000, -109701/10240000, 51013/10240000, -17/12800, 1/6400 ], P![]> *];
models[[Integers()|1,10,122,305]] := [* <1, P![ 0, 625, -46875/16, 25625/8, -16875/16 ], P![]> *];
models[[Integers()|1,2,61,122]] := [* <3, P![ 0, 1953125/16, -352734375/256, 1572265625/256, -14188671875/1024, 8741796875/512, -12055859375/1024, 137890625/32, -10546875/16 ], P![]> *];
models[[Integers()|1,5,61,305]] := [* <3, P![ -15625, 703125/4, -202796875/256, 481078125/256, -2680484375/1024, 1141578125/512, -1176734375/1024, 84765625/256, -10546875/256 ], P![]> *];
models[[Integers()|1,2,5,10]] := [* <5, P![ -48828125/16, 3486328125/64, -1739990234375/4096, 1957490234375/1024, -45324326171875/8192, 89151728515625/8192, -982382099609375/65536, 480123193359375/32768, -664967060546875/65536, 79986259765625/16384, -25513193359375/16384, 152099609375/512, -6591796875/256 ], P![]> *];
