// CM points as fixed points of quadratic elements of the order, and the cross-ratio used to compare
// a computed Hauptmodul with a published one.  Load after vw3.m and vwal.m.
//
// A CM point of discriminant d on X_0^D(N) is the fixed point in the upper half plane of an element
// mu of O generating a quadratic order of discriminant d, i.e. with reduced trace t and reduced
// norm n satisfying t^2 - 4n = d.  Taking the fixed point of mu directly avoids Magma's CMPoints,
// which is broken at these discriminants (RingClassGroup coercion error inside
// package/Geometry/GrpPSL2/GrpPSL2Shim/reciprocity.m line 29).
//
// A computed Hauptmodul is only defined up to a Mobius change of coordinate, so it is compared with
// a published one through CROSS-RATIOS against three fixed reference points.  Nothing is fitted:
// three points set the coordinate, and every further value is then a prediction.

// An element of O of reduced trace tr and reduced norm nr.  The reduced trace is a linear form on
// the coefficients, so one coefficient is solved for and only the other three are searched.
function ElementTraceNorm(O, tr, nr, B)
    bas := Basis(O);
    ts := [Integers() | Trace(b) : b in bas];
    js := [j : j in [1..4] | ts[j] ne 0];
    if IsEmpty(js) then return false, _; end if;
    j := js[1];
    rest := [i : i in [1..4] | i ne j];
    for c1 in [-B..B] do for c2 in [-B..B] do for c3 in [-B..B] do
        cs := [c1, c2, c3];
        num := tr - &+[ cs[i]*ts[rest[i]] : i in [1..3] ];
        if num mod ts[j] ne 0 then continue; end if;
        co := [Integers() | 0, 0, 0, 0];
        for i in [1..3] do co[rest[i]] := cs[i]; end for;
        co[j] := num div ts[j];
        mu := &+[ co[i]*bas[i] : i in [1..4] ];
        if Norm(mu) eq nr then return true, mu; end if;
    end for; end for; end for;
    return false, _;
end function;

// (trace, norm) pairs with trace^2 - 4*norm = d and 0 < norm <= nmax, smallest norm first: those
// elements have the smallest coefficients and are found in the smallest box.
function TraceNormPairs(d, nmax)
    out := [];
    for t in [0..Isqrt(4*nmax - d)] do
        if (t^2 - d) mod 4 ne 0 then continue; end if;
        n := (t^2 - d) div 4;
        if n gt 0 and n le nmax then Append(~out, <t, n>); end if;
    end for;
    return out;
end function;

// A CM point of discriminant d, as a disc point reduced into the fundamental domain.
//
// Magma's Embed(R, O) returns an OPTIMAL embedding of the quadratic order R into O, so the fixed
// point really is a CM point of discriminant disc(R) and not of some larger order containing it --
// which a naive search for an element of given trace and norm does not guarantee.  It also finds
// elements a box search misses: at D = 10 the generator of discriminant -40 has coefficients
// outside a box of 6, and a brute-force search over that box returns nothing at all.
function CMDiscPoint(O, d, DD, mrep, p, CC, G)
    R := QuadraticOrder(BinaryQuadraticForms(d));
    ok := true; mu := O ! 1;
    try mu := Embed(R, O); catch e ok := false; end try;
    if not ok then return false, _, _, _; end if;
    if Trace(mu)^2 - 4*Norm(mu) ne d then return false, _, _, _; end if;
    W := DiscMatrix(mrep, mu, CC ! p, CC);
    a := W[1][1]; b := W[1][2]; c := W[2][1]; e := W[2][2];
    if Abs(c) lt 10^-30 then return false, _, _, _; end if;
    rad := Sqrt((e-a)^2 + 4*b*c);
    for rt in [ (-(e-a)+rad)/(2*c), (-(e-a)-rad)/(2*c) ] do
        if Abs(rt) lt 0.9999 then
            red := ShimuraReduceUnit(G ! (O ! 1), Generators(G), DD : z0 := DD ! rt);
            Wd := DiscMatrix(mrep, Quaternion(red[1]), CC ! p, CC);
            return true, ActOn(Wd, rt), mu, <Trace(mu), Norm(mu)>;
        end if;
    end for;
    return false, _, _, _;
end function;

// cross-ratio (z1, z2 ; z3, z4), with Infinity() allowed in any one argument.
function CrossRatio(z1, z2, z3, z4)
    if Type(z1) eq Infty then return (z2-z4)/(z2-z3); end if;
    if Type(z2) eq Infty then return (z1-z3)/(z1-z4); end if;
    if Type(z3) eq Infty then return (z2-z4)/(z1-z4); end if;
    if Type(z4) eq Infty then return (z1-z3)/(z2-z3); end if;
    return ((z1-z3)*(z2-z4))/((z1-z4)*(z2-z3));
end function;
