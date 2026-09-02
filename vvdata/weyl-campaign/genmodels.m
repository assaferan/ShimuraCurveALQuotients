// Model-generation harness, recreated from HANDOFF.md (the original lived at /tmp/genmodels.m
// and was lost to a nightly purge).  Two deliberate changes from the HANDOFF version:
//   * no `sage -sh` wrapper -- the polymake backend is dead and replaced by nmzsolve.py,
//     which needs only python3 + NORMALIZ_BIN;
//   * the output directory is a parameter (OUTDIR), defaulting to the scratchpad, so a
//     probe run does not mutate data/models/ in the repo.
//
//   NORMALIZ_BIN=... magma -b D_s:=58 N_s:=3 OUTDIR:=/path genmodels.m
SetQuitOnError(true);
SetColumns(0);
AttachSpec("ShimuraQuotients.spec");
// Progress ON by default.  The [n/6] stage markers in EquationsCovers.m are vprintf level 1,
// so without this a multi-hour run logs nothing until it ends and "is it progressing or stuck?"
// is unanswerable.  That cost visibility four separate times on 2026-08-31.  VERB:=0 to mute.
verb := 1;
if assigned VERB then verb := StringToInteger(VERB); end if;
SetVerbose("ShimuraQuotients", verb);
Px<x> := PolynomialRing(Rationals());
D := StringToInteger(D_s); N := StringToInteger(N_s);
outdir := "data/models";
if assigned OUTDIR then outdir := OUTDIR; end if;
// Skip known-deferred vx bases (memory vx-laurent-n0-circular): huge Zero-side space then
// crash on the AbsEltseq "vx ge 0" assert. Fail fast.
vx_skip := {<95,1>, <115,1>, <123,1>, <129,1>};
if <D,N> in vx_skip then WriteStderr(Sprintf("SKIP vx base %o_%o\n", D, N)); quit; end if;
curves := GetHyperellipticCandidates();
Xstar := rep{X : X in curves | X`D eq D and X`N eq N and IsStarCurve(X)};
t0 := Realtime();
// INTSOL=1 turns on IntegralSolution for an UNCAPPED run.  gencapped.m hardwires it together
// with the genus cap, so this is the way to vary integrality alone and keep every cover --
// which is what separates "integrality rescued this base" from "the genus cap did".
intsol := GetEnv("INTSOL") ne "";
WriteStderr(Sprintf("IntegralSolution = %o\n", intsol));
covers, ws := AllEquationsAboveCovers(Xstar, curves : Prec := 100, IntegralSolution := intsol);
WriteStderr(Sprintf("computed AllEquationsAboveCovers in %os\n", Realtime()-t0));
agg := AssociativeArray();
for label in Keys(covers) do
  X := curves[label]; Wkey := Sort([w : w in X`W]);
  if not IsDefined(agg, Wkey) then agg[Wkey] := [* *]; end if;
  for base in Keys(covers[label]) do Append(~agg[Wkey], <X`g, covers[label][base]>); end for;
end for;
keystr := func<W | "[Integers()|" cat (#W eq 0 select "" else &cat[IntegerToString(w) cat (i lt #W select "," else "") : i->w in W]) cat "]">;
lines := [
  Sprintf("// Subhyperelliptic cover models for X_0(%o,%o)*", D, N),
  "// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).",
  "P<x> := PolynomialRing(Rationals());",
  "models := AssociativeArray();"
];
ncov := 0;
for Wkey in Keys(agg) do
  entries := [];
  for m in agg[Wkey] do
    g := m[1]; C := m[2];
    if Type(C) eq CrvHyp then
      f, h := HyperellipticPolynomials(C);
      Append(~entries, Sprintf("<%o, P!%o, P!%o>", g, Coefficients(f), Coefficients(h)));
    else
      Append(~entries, Sprintf("<%o, \"CRV\", %m>", g, [Sprint(p) : p in DefiningPolynomials(C)]));
    end if;
  end for;
  Append(~lines, Sprintf("models[%o] := [* %o *];", keystr(Wkey), Join(entries, ", ")));
  ncov +:= 1;
end for;
fname := Sprintf("%o/models_%o_%o.m", outdir, D, N);
Write(fname, Join(lines, "\n") : Overwrite);
WriteStderr(Sprintf("WROTE %o (%o cover-keys)\n", fname, ncov));
printf "GENMODELS_DONE %o %o keys %o\n", D, N, ncov;
quit;
