// Subhyperelliptic cover models for X_0(14,43)*
//
// PROVENANCE. Produced 2026-09-06 on lovelace; AllEquationsAboveCovers took 151189 s (42 h).
// 4 cover-keys, all populated (genus 1, 3, 2, 6).
//
// ⚠ PRODUCED WITH INTSOL=1 -- this file is NOT known to regenerate under default settings.
//     INTSOL=1 NORMALIZ_BIN=... magma -b D_s:=14 N_s:=43 OUTDIR:=... genmodels.m < /dev/null
// ⚠ CONFIDENCE IN THAT FLAG: this session captured the launch wrapper in `ps`, which showed
// 10_61 and 14_43 started from one command line with `INTSOL=1` prefixed to the 14_43 half only.
// The run LOG does not independently record it (lovelace's genmodels.m predates the line that
// prints `IntegralSolution = ...`), and the wrapper has since exited. So the flag rests on that
// one observation -- treat it as the best available record, not as a log-confirmed fact, and if
// this file is ever regenerated, try INTSOL=1 first.
//
// WHY THIS ONE MATTERS. 14_43 is from the OBSTRUCTED class, the second base out of it (after
// 10_61). It is also one of only two bases that triage found RUNNABLE among 122 never-started
// ones.
//
// ⚠ NOT a Guo-Yang base -- X_0^14(43) has no published equation, so there is NO external oracle.
// ModelChecks alone validates it; its strongest component is the Eichler-Selberg point count,
// which is independent of the Borcherds/Schofer path that produced the model. Weaker evidence
// than the published-equation bases carry, and it should be quoted that way.
//
// ⚠ CODE VERSION / PROCESS HAZARD. Started 2026-09-04 from lovelace's clone, which was `git
// pull`-ed to main ON 2026-09-06 WHILE THIS JOB WAS STILL RUNNING FROM IT. AttachSpec loads
// packages on demand, so the output cannot be pinned to a single commit. No harm is evident
// (ModelChecks passes), but see data/models/PROVENANCE.md -- launch long runs from a COPIED tree.
//
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,7,86,602]] := [* <1, P![ -176, 1440, -5040, 6912 ], P![]> *];
models[[Integers()|1,2,43,86]] := [* <3, P![ -23068672, 119537664, -258195456, 84492288, -1540435968, 6372089856, -6064533504, 1576599552, -13759414272 ], P![]> *];
models[[Integers()|1,2,301,602]] := [* <2, P![ 16384, -172032, 746496, -1714176, 1741824, -995328, 2985984 ], P![]> *];
models[[Integers()|1,2,7,14]] := [* <6, P![ -377957122048, 5927054868480, -42015249137664, 174580354252800, -477605916573696, 1106244676878336, -3127921531158528, 9227540203831296, -19484926701207552, 26824253311549440, -34479303441776640, 51395375017230336, -43644311694213120, 18402831325200384, -41085390865563648 ], P![]> *];
