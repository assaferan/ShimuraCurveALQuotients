cached_orders := NewStore();
cached_traces := NewStore();
class_nos := NewStore();
point_counts := NewStore();

// Collect mode: a "dry run" of a trace computation in which the leaf class-number function
// records the discriminants it would need (so they can be fetched in one batched, sorted
// streaming pass) instead of computing them.  While collecting, the memoization caches are
// bypassed -- GetCache always misses and SetCache is a no-op -- so the dry run fully traverses
// the summation (calling H to record discriminants) and leaves no dummy values cached.
collecting := NewStore();

procedure StartCollecting()
  StoreSet(collecting, "on", true);
  StoreSet(collecting, "discs", AssociativeArray());
end procedure;

procedure StopCollecting()
  StoreSet(collecting, "on", false);
end procedure;

function IsCollecting()
  bool, v := StoreIsDefined(collecting, "on");
  return bool and v;
end function;

procedure CollectDisc(D)
  bool, s := StoreIsDefined(collecting, "discs");
  if not bool then s := AssociativeArray(); end if;
  s[D] := true;                       // associative-array key gives free dedup
  StoreSet(collecting, "discs", s);
end procedure;

function GetCollectedDiscs()
  bool, s := StoreIsDefined(collecting, "discs");
  if not bool then return []; end if;
  return Setseq(Keys(s));
end function;

intrinsic CacheClear(name)
{Clear the internal cache for cached_orders}
  // We need to save and restore the id, otherwise horrific things might
  // happen
  StoreClear(name);
  StoreSet(name, "cache", AssociativeArray());
end intrinsic;

procedure SetCache(k,v, name)
  if IsCollecting() then return; end if;     // don't cache dummy dry-run values
  bool, cache := StoreIsDefined(name, "cache");
  if not bool then
    cache := AssociativeArray();
  end if;
  cache[k] := v;
  StoreSet(name, "cache", cache);
end procedure;

function GetCache(k, name)
  if IsCollecting() then return false, _; end if;   // force full traversal while collecting
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
