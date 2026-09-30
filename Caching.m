cached_orders := NewStore();
cached_traces := NewStore();
class_nos := NewStore();
point_counts := NewStore();
// data/hypg<g>q<p>.txt tables loaded by CheckWeilPolynomialAtPrime, keyed by <g, p>
weil_tables := NewStore();
// SpecialFiberIsomorphism source certificates, keyed by <D, N, W, g, p>
sfi_certificates := NewStore();
// DualGraphData fibers of X_0^D(N) at p | D, keyed by <D, N, p, cache dir>
dual_graphs := NewStore();

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
  if not bool then s := AssociativeArray(); else StoreRemove(collecting, "discs"); end if;   // in place, see SetCache
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

// ⚠ Drop the store's reference BEFORE mutating.  Magma copies a value on write when more than
// one reference to it exists; with the store still holding the array, `cache[k] := v` copied
// the whole associative array on EVERY insert -- quadratic in cache size.  Measured on a Mac
// (2026-09-30, vvdata/weyl-campaign/weil-cost-2026-09-30/cache_bench{,2}.m on the campaign
// branch): 20k inserts 4.3 s, 80k inserts 68 s; 1.28M inserts 2.4 s once the store's reference
// is removed first.  The class-number cache grows to hundreds of thousands of entries across a
// Weil-stage run.  Between the StoreRemove and the StoreSet the array lives only in the local
// `cache`: an error raised there loses the cache (a recomputation cost, never a wrong value).
procedure SetCache(k,v, name)
  bool, cache := StoreIsDefined(name, "cache");
  if not bool then
    cache := AssociativeArray();
  else
    StoreRemove(name, "cache");
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
