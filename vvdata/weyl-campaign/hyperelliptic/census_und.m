// For each undecided curve: does its base have a model file, and do the stored double covers above
// it generate its function field (so that a fibre product gives a model)?
AttachSpec("ShimuraQuotients.spec"); SetColumns(0);
P<x> := PolynomialRing(Rationals());
curves := eval Read("data/curves_after_UpdateCurves8.dat");
und := [X : X in curves | not assigned X`IsSubhyp];
nofile := 0; buildable := []; notbuildable := []; stored := 0;
files := AssociativeArray();
for X in und do
    fn := Sprintf("data/models/models_%o_%o.m", X`D, X`N);
    if not FileExists(fn) then nofile +:= 1; continue; end if;
    if not IsDefined(files, fn) then files[fn] := eval (Read(fn) cat "\nreturn models;"); end if;
    models := files[fn];
    k := Sort(SetToSequence(X`W));
    if IsDefined(models, k) and #models[k] gt 0 then stored +:= 1; Append(~buildable, <X`D, X`N, k, X`g, "stored">); continue; end if;
    hyp := {Set(kk) : kk in Keys(models) | #models[kk] gt 0 and Type(models[kk][1][2]) ne MonStgElt};
    ups := AtkinLehnerDoubleCoversOver(X`W, X`D, X`N);
    avail := {U : U in ups | U in hyp};
    ok, gens := FibreProductGenerators(X`W, avail, X`D, X`N);
    if ok then Append(~buildable, <X`D, X`N, k, X`g, "fibre product">); else Append(~notbuildable, <X`D, X`N, k, X`g, #avail, #ups>); end if;
end for;
printf "undecided %o: no model file for the base %o; entry already stored %o; buildable as fibre product %o; not buildable %o\n",
    #und, nofile, stored, #[b : b in buildable | b[5] eq "fibre product"], #notbuildable;
printf "buildable: %o\n", buildable;
printf "not buildable (D, N, W, g, stored double covers above W, index-2 overgroups): %o\n", notbuildable;
quit;
