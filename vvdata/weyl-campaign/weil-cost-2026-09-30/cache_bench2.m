// Candidate fix: drop the store's reference BEFORE mutating the local copy, so the associative
// array has refcount 1 and the insert happens in place (no copy-on-write), then put it back.
AttachSpec("ShimuraQuotients.spec");
st := NewStore();
procedure SetCacheFast(k, v, name)
    bool, cache := StoreIsDefined(name, "cache");
    if not bool then cache := AssociativeArray(); else StoreRemove(name, "cache"); end if;
    cache[k] := v;
    StoreSet(name, "cache", cache);
end procedure;
function GetCacheFast(k, name)
    bool, cache := StoreIsDefined(name, "cache");
    if not bool then return false, _; end if;
    return IsDefined(cache, k);
end function;
for Nins in [20000, 80000, 320000, 1280000] do
    StoreClear(st);
    t0 := Cputime();
    for i in [1..Nins] do SetCacheFast(-4*i, i, st); end for;
    t_ins := Cputime(t0);
    t0 := Cputime();
    hits := 0;
    for i in [1..Nins] do b, v := GetCacheFast(-4*i, st); if b then hits +:= 1; end if; end for;
    t_get := Cputime(t0);
    printf "N=%o  SetCacheFast %o s   GetCache %o s   hits %o\n", Nins, RealField(4)!t_ins, RealField(4)!t_get, hits;
end for;
exit;
