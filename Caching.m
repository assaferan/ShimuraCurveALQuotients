cached_orders := NewStore();
cached_traces := NewStore();
class_nos := NewStore();

intrinsic CacheClear(name)
{Clear the internal cache for cached_orders}
  // We need to save and restore the id, otherwise horrific things might
  // happen
  StoreClear(name);
  StoreSet(name, "cache", AssociativeArray());
end intrinsic;

procedure SetCache(k,v, name)
  bool, cache := StoreIsDefined(name, "cache");
  if not bool then
    cache := AssociativeArray();
  end if;
  cache[k] := v;
  StoreSet(name, "cache", cache);
end procedure;

function GetCache(k, name)
  bool, cache := StoreIsDefined(name, "cache");
  if not bool then
    cache := AssociativeArray();
    return false, _;
  end if;
  return IsDefined(cache, k);
end function;


// intrinsic BinomialCached(n::RngIntElt, k::RngIntElt) -> RngIntElt
// {The binomial coefficient n choose r}
//     b, v := GetCache(n, k);
//     if not b then
//         "caching", <n,k>;
//         v := Binomial(n, k);
//         SetCache(n ,k, v);
//     else
//         "already computed it!";
//     end if;
//     return v;
// end intrinsic;
