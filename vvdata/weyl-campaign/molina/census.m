// Which hyperelliptic Atkin-Lehner quotients at level 1 with odd discriminant have no model?
// Models are scattered across branches, so the search is over every data/models directory found
// in the checkout, unioned.
//   magma -b modeldir:=... logfile:=... census.m < /dev/null

AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
L := Open(logfile, "w");
procedure Say(s) fprintf L, "%o\n", s; Flush(L); end procedure;

// every models_D_1.m anywhere under modeldir, keyed by <D, sorted W>
have := {};
files := Split(Pipe(Sprintf("find %o -name 'models_*_1.m'", modeldir), ""), "\n");
nf := 0;
for fn in files do
    if fn eq "" then continue; end if;
    nf +:= 1;
    bits := Split(Split(fn, "/")[#Split(fn, "/")], "_.");
    DD := StringToInteger(bits[2]);
    try
        models := eval (Read(fn) cat "\nreturn models;");
    catch err
        continue;
    end try;
    for k in Keys(models) do
        if #models[k] eq 0 then continue; end if;
        if Type(models[k][1][2]) eq MonStgElt then continue; end if;   // fibre-product placeholder
        Include(~have, <DD, Sort(k)>);
    end for;
end for;
Say(Sprintf("read %o model files, %o <D, W> keys with a y^2 = f model", nf, #have));

curves := GetHyperellipticCandidates();
miss := [];
for X in curves do
    if X`N ne 1 or IsEven(X`D) or X`D eq 1 then continue; end if;
    if not assigned X`g or X`g lt 2 then continue; end if;
    if not (assigned X`IsSubhyp and X`IsSubhyp and assigned X`IsHyp and X`IsHyp) then continue; end if;
    k := Sort(SetToSequence(X`W));
    if <X`D, k> in have then continue; end if;
    Append(~miss, <X`D, k, X`g>);
end for;
miss := Sort(miss, func<a,b | a[3] ne b[3] select a[3]-b[3] else (a[1] ne b[1] select a[1]-b[1] else 0)>);
Say(Sprintf("MISSING: %o hyperelliptic quotients with odd D, N = 1, genus >= 2 and no model", #miss));
for t in miss do
    Say(Sprintf("  D = %o  W = %o  genus %o  (#W = %o)", t[1], t[2], t[3], #t[2]));
end for;
bygenus := AssociativeArray();
for t in miss do
    if IsDefined(bygenus, t[3]) then bygenus[t[3]] +:= 1; else bygenus[t[3]] := 1; end if;
end for;
Say(Sprintf("by genus: %o", [<k, bygenus[k]> : k in Sort(Setseq(Keys(bygenus)))]));
Say("CENSUS DONE");
quit;
