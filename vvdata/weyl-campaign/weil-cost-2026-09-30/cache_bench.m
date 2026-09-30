// Is the SetCache pattern (StoreIsDefined -> local copy -> insert -> StoreSet) quadratic?
// Compare N inserts through Caching.m's SetCache against a plain associative array.
AttachSpec("ShimuraQuotients.spec");
_ := ClassNumberLU(-4);
import "Caching.m" : SetCache, GetCache, class_nos;
for Nins in [20000, 40000, 80000] do
    CacheClear(class_nos);
    t0 := Cputime();
    for i in [1..Nins] do SetCache(-4*i, i, class_nos); end for;
    t_store := Cputime(t0);
    A := AssociativeArray();
    t0 := Cputime();
    for i in [1..Nins] do A[-4*i] := i; end for;
    t_plain := Cputime(t0);
    printf "N=%o  SetCache %o s   plain assoc %o s   ratio %o\n", Nins, RealField(4)!t_store, RealField(4)!t_plain,
        t_plain gt 0 select RealField(4)!(t_store/t_plain) else 0;
end for;
exit;
