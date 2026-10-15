// Atkin-Lehner layer on top of vw3.m (load vw3.m first).
//
// The published method expands forms on X_0^D(N) itself.  Every model this project stores is of a
// quotient X_0^D(N)/W, so the group has to be extended by the Atkin-Lehner elements and the space
// of weight-k forms cut down to the W-invariant part.
//
// The Atkin-Lehner element for a Hall divisor m of DN is an element mu of the Eichler order O with
// reduced norm m that normalises O.  Such a mu generates the two-sided ideal of reduced norm m, and
// conversely every element of that ideal of reduced norm exactly m is a generator and so normalises
// O.  Searching the ideal's lattice therefore finds one quickly, where searching O does not: a box
// of 3 over the ideal basis succeeds at every base tried, while a box of 4 over O misses several.
// The reduced norm MUST come out positive; a generator of negative norm maps the upper half plane
// to the lower one and does not act on the curve (Magma's IsPrincipal returns such a generator at
// D = 6, m = 3, so taking the first generator it offers is wrong).
//
// Magma's Normalizer is implemented for definite orders only, and these orders are indefinite, so
// it cannot be used here.

function ALNormalises(mu, O)
    for b in Basis(O) do if not (mu*b*mu^-1 in O) then return false; end if; end for;
    return true;
end function;

// mu in O with reduced norm exactly m (positive) normalising O, found in the two-sided ideal.
function ALElement(O, m, B)
    P := 1*O;
    for p in PrimeDivisors(m) do P := P * PrimeIdeal(O, p); end for;
    bas := Basis(P);
    for c1 in [-B..B] do for c2 in [-B..B] do for c3 in [-B..B] do for c4 in [-B..B] do
        mu := c1*bas[1] + c2*bas[2] + c3*bas[3] + c4*bas[4];
        if mu eq 0 or Norm(mu) ne m then continue; end if;
        if ALNormalises(mu, O) then return true, mu; end if;
    end for; end for; end for; end for;
    return false, _;
end function;

// Atkin-Lehner product of Hall divisors of DN.
function ALMul(u, v, DN)
    g := GCD(u, v);
    return (u div g) * (v div g);
end function;

// A minimal generating set of the subgroup W (a set of Hall divisors containing 1).
function ALGens(W, DN)
    span := {Integers() | 1}; gens := [];
    for m in Sort(SetToSequence(W)) do
        if m in span then continue; end if;
        Append(~gens, m);
        span := {ALMul(s, m, DN) : s in span} join span;
    end for;
    return gens;
end function;

// Mobius action of a 2x2 matrix on a disc point.
function ActOn(W, w) return (W[1][1]*w + W[1][2]) / (W[2][1]*w + W[2][2]); end function;

// Samples carrying the Atkin-Lehner relations, in the same <source, target, matrix> format that
// ReducedSamples produces, so AutSystem consumes both without change.  For a sample point wc on the
// circle of radius rad we apply mu, then reduce the image back into the fundamental domain by an
// element delta of Gamma; the relation imposed is the one for the composite delta*mu, which lies in
// the extended group.
//
// ⚠ ShimuraReduceUnit does not always succeed here.  Applying mu throws the point very close to the
// boundary (modulus up to 0.997 at D = 6), and the reduction then returns a point still outside the
// domain -- 20 to 30 percent of them, with modulus up to 0.91 against a domain radius of 0.54.  The
// truncated series is meaningless there (0.91^50 is about 1 percent), and keeping those rows
// destroys the linear system.  They are dropped.  On the rows that do reduce, the relation holds to
// the same accuracy as the Gamma relation itself (3e-10 against 3e-11 at D = 6, weight 4), so the
// relation is right and only the reduction is at fault.
function ALSamples(G, O, DD, mrep, p, rad, Q, CC, mus)
    out := []; gens := Generators(G);
    pi := Pi(RealField(Precision(CC)));
    tol := rad + 1e-10;
    nbad := 0;
    for mu in mus do
        Wmu := DiscMatrix(mrep, mu, CC ! p, CC);
        for m in [1..Q] do
            wc := CC ! (rad * Exp(2*pi*CC.1*m/Q));
            wmu := ActOn(Wmu, wc);
            if Abs(wmu) ge 1 then continue; end if;
            red := ShimuraReduceUnit(G ! (O ! 1), gens, DD : z0 := DD ! wmu);
            Wd := DiscMatrix(mrep, Quaternion(red[1]), CC ! p, CC);
            Wt := Wd * Wmu;
            wt := ActOn(Wt, wc);
            if Abs(wt) gt tol then nbad +:= 1; continue; end if;
            Append(~out, <wc, wt, Wt>);
        end for;
    end for;
    return out, nbad;
end function;

// The signature of the quotient X = X_0^D(N)/W, from the repository's own elliptic-point counts:
// genus from the genus formula, and the elliptic orders with multiplicity.
function QuotientSignature(X)
    g := GenusShimuraCurveQuotient(X`D, X`N, X`W);
    ell := NumberOfEllipticPointsByCMOrder(X);
    es := [];
    for q in Sort(SetToSequence(Keys(ell))) do
        n := &+([0] cat [ell[q][d] : d in Keys(ell[q])]);
        es cat:= [q : i in [1..n]];
    end for;
    return g, es;
end function;

// dim M_k for a cocompact Fuchsian group of signature (g; e_1..e_r): the Riemann-Roch count for
// k >= 4 even, and the genus itself at k = 2 (weight-2 forms are the holomorphic differentials).
function DimMkSig(k, g, es)
    if k eq 0 then return 1; end if;
    if k eq 2 then return g; end if;
    return (k-1)*(g-1) + &+([0] cat [Floor((k/2)*(1 - 1/e)) : e in es]);
end function;

// Number of terms needed for a truncation floor of about 10^-target at domain radius rho, and the
// absolute level below which a singular value counts as vanishing.  Both must follow rho: at
// D = 6 rho is 0.54 and 50 terms suffice, but at D = 10 rho is 0.82 and 50 terms leave an error of
// 1e-5, so a fixed tolerance of 1e-7 reports every space as empty.
function TermsFor(rho, target)
    return Ceiling(target * Log(10) / (-Log(rho)));
end function;
function ZeroLevel(rho, nterms)
    return Max(1e-10, 1e6 * rho^nterms);
end function;

// Dimension as the position of the largest gap in the ascending singular values, counting only
// values that are small enough to be plausibly zero.  Returns the dimension and the gap ratio.
function GapDim(sv, tolabs)
    nz := #[s : s in sv | s lt tolabs];
    if nz eq 0 then return 0, 0.0; end if;
    best := 0.0; bi := 0;
    for i in [1..Min(nz, #sv-1)] do
        if sv[i] gt 0 and sv[i+1]/sv[i] gt best then best := sv[i+1]/sv[i]; bi := i; end if;
    end for;
    return bi, best;
end function;
