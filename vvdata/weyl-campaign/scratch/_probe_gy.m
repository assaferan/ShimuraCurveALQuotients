// Which Guo-Yang bases do we have models for, and which are missing?
// The GY target list = every base with a GuoYang CM-table offline test, plus the equation tables.
gy := Split(Pipe("ls tests/_offline/GuoYang_*.m 2>/dev/null | sed 's|.*GuoYang_||; s|\\.m||'", ""), "\n");
gy := [g : g in gy | #g gt 0];
printf "GuoYang offline CM-table tests: %o\n", #gy;
have := []; miss := [];
for g in gy do
    fn := "data/models/models_" cat g cat ".m";
    if not FileExists(fn) then Append(~miss, g cat "(no file)"); continue; end if;
    models := eval (Read(fn) cat "\nreturn models;");
    pop := [k : k in Keys(models) | #models[k] gt 0];
    tot := #Keys(models);
    k1 := [Integers()|1];
    hasfull := IsDefined(models, k1) and #models[k1] gt 0;
    if #pop eq 0 then
        Append(~miss, g cat "(all keys empty)");
    elif #pop lt tot then
        Append(~miss, Sprintf("%o(%o/%o keys populated%o)", g, #pop, tot, hasfull select "" else ", NO W=[1]"));
    else
        Append(~have, g);
    end if;
end for;
printf "\nFULLY POPULATED (%o): %o\n", #have, have;
printf "\nINCOMPLETE OR MISSING (%o):\n", #miss;
for m in miss do printf "   %o\n", m; end for;
