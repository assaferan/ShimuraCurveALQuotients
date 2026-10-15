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
//
// The cube in the coordinates of a basis of that ideal is the wrong shape for the same reason as
// in the Hecke coset search: the reduced norm form is indefinite, so its level set Nrd = m is a
// hyperboloid and a cube cuts it along directions that have nothing to do with the lattice.  A box
// of 3 finds nothing at all at D = 35 for m = 7.  Enumerating short vectors of the positive
// definite majorant of vw3.m instead, and raising the bound until an element turns up, finds one
// at every base tried; the second argument is kept so existing callers still work, and is used
// only as a starting bound.
// ⚠ A SCALAR IS NOT AN ATKIN-LEHNER ELEMENT, and at non-squarefree level the search can return
// one.  For m a square, say m = 9 at level 9, the element 3 of the order has reduced norm 9 and
// normalises the order, but it acts trivially on the upper half plane, so taking it as w_9 builds
// a group SMALLER than the normaliser and the invariant space comes out too big -- at X_0^10(9)
// the weight-4 space on the star had dimension 4 where the signature predicts 3, and the pencil
// never resolved.  Reject anything with mu/d still in the order.  The search also runs over the
// order itself rather than over a product of PrimeIdeal(O, p), which is the wrong ideal when a
// prime divides the level more than once.
function ALElement(O, m, B : cmax := 65536)
    A := Algebra(O);
    mrep := FuchsianMatrixRepresentation(A : Precision := 60);
    bas := Basis(O);
    L := LatticeWithGram(FrobeniusGram(bas, mrep, A));
    ds := [ d : d in [2..Isqrt(m)] | m mod d^2 eq 0 ];
    c := 8;
    while c le cmax do
        for v in ShortVectors(L, SCL*c*m) do
            mu := &+[ (Integers() ! v[1][i]) * bas[i] : i in [1..4] ];
            if Norm(mu) ne m then continue; end if;
            if exists{ d : d in ds | (A ! mu)/d in O } then continue; end if;
            if ALNormalises(mu, O) then return true, mu; end if;
        end for;
        c *:= 4;
    end while;
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
// The reduction back into the domain descends over the SIDE PAIRING, for the same reason as in
// vw3.m: ShimuraReduceUnit given only Generators(G) returns the identity on the points it cannot
// handle, and here it failed on 20 to 30 percent of the Atkin-Lehner images even at D = 6, leaving
// them at modulus up to 0.91 against a domain radius of 0.54, where a truncated series is
// worthless.  Applying mu throws a sample point very close to the boundary, so this is the harder
// of the two reductions and the one where the defect showed first.  Rows whose target still sits
// outside the sampling radius are dropped and counted, so a caller can see whether any were.
function ALSamples(G, O, DD, mrep, p, rad, Q, CC, mus)
    out := [];
    Ws := SidePairingMatrices(G, DD, mrep, p, CC);
    pi := Pi(RealField(Precision(CC)));
    tol := rad + 1e-10;
    nbad := 0;
    for mu in mus do
        Wmu := DiscMatrix(mrep, mu, CC ! p, CC);
        for m in [1..Q] do
            wc := CC ! (rad * Exp(2*pi*CC.1*m/Q));
            wmu := ActOn(Wmu, wc);
            if Abs(wmu) ge 1 then continue; end if;
            Wd, _ := ReduceByPairing(Ws, wmu, CC);
            Wt := Wd * Wmu;
            wt := ActOn(Wt, wc);
            if Abs(wt) gt tol then nbad +:= 1; continue; end if;
            Append(~out, <wc, wt, Wt>);
        end for;
    end for;
    return out, nbad;
end function;

// ---- expanding on the NORMALISER's domain instead of Gamma's ----
//
// The term count is NT = T ln10 / (-ln rho), with rho the circumradius of the domain the points
// are reduced into, and the decomposition costs NT^3.  Every model here is of a quotient by an
// Atkin-Lehner subgroup, so the forms wanted are invariant under the larger group generated by
// Gamma and the Atkin-Lehner elements; its domain is 1/#W* of the area and has a visibly smaller
// circumradius.  At X_0^35(1) rho falls from 0.879 to 0.774 and the term count from 357 to 180, a
// factor of seven in the decomposition; the gain grows as rho approaches 1, which is exactly the
// regime where the cost is otherwise prohibitive.
//
// Reducing over the larger group is a two-level descent: reduce by the side pairing of Gamma,
// then try each Atkin-Lehner translate and reduce by Gamma again, keeping whichever decreases the
// modulus most, and repeat.  The Atkin-Lehner steps are recorded as a Hall divisor, so each
// sample comes with the coset of the larger group it landed in; that label is what lets ONE set
// of samples serve the invariant space, the star space and the anti-invariant space at once,
// since a relation for an element labelled m holds with a sign that depends only on m.
function ReduceOverStar(Ws, Als, w, CC : maxit := 40)
    eps := (RealField(Precision(CC)) ! 10)^(-Precision(CC) div 3);
    Acc, cur := ReduceByPairing(Ws, w, CC);
    lab := 1;
    for it in [1..maxit] do
        best := Abs(cur); bM := Acc; bw := cur; bl := lab; moved := false;
        for am in Als do
            m := am[1]; M := am[2];
            den := M[2][1]*cur + M[2][2];
            if Abs(den) le eps then continue; end if;
            nw := (M[1][1]*cur + M[1][2])/den;
            if Abs(nw) ge 1 then continue; end if;
            A2, nw2 := ReduceByPairing(Ws, nw, CC);
            if Abs(nw2) lt best - eps then
                best := Abs(nw2); bM := A2*M*Acc; bw := nw2;
                bl := (lab*m) div GCD(lab, m)^2; moved := true;
            end if;
        end for;
        if not moved then break; end if;
        Acc := bM; cur := bw; lab := bl;
    end for;
    return Acc, cur, lab;
end function;

// Samples on the circle of radius rad, reduced over the larger group, as <w, w', matrix, label>.
// A point that reduces to itself carries no relation and is dropped.
function StarSamples(Ws, Als, rad, Q, CC)
    RR := RealField(Precision(CC));
    eps := (RR ! 10)^(-Precision(CC) div 3);
    out := [];
    for m in [1..Q] do
        wc := CC ! (rad * Exp(2*Pi(RR)*CC.1*m/Q));
        Acc, wp, lab := ReduceOverStar(Ws, Als, wc, CC);
        if Abs(wp - wc) le eps then continue; end if;
        Append(~out, <wc, wp, Acc, lab>);
    end for;
    return out;
end function;

// The same, but each point is first moved by the Atkin-Lehner element mu of Hall divisor m and
// then reduced; the relation obtained is the one for an element of the coset labelled m times
// whatever the reduction used.  This is what supplies rows in the cosets the plain sampling
// happens not to reach.
function StarSamplesVia(Ws, Als, Wmu, m0, rad, Q, CC)
    RR := RealField(Precision(CC));
    out := [];
    for i in [1..Q] do
        wc := CC ! (rad * Exp(2*Pi(RR)*CC.1*i/Q));
        den := Wmu[2][1]*wc + Wmu[2][2];
        nw := (Wmu[1][1]*wc + Wmu[1][2])/den;
        if Abs(nw) ge 1 then continue; end if;
        Acc, wp, lab := ReduceOverStar(Ws, Als, nw, CC);
        if Abs(wp) gt rad then continue; end if;
        Append(~out, <wc, wp, Acc*Wmu, (lab*m0) div GCD(lab, m0)^2>);
    end for;
    return out;
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
// Measured separation, D = 6 and D = 10 with the disc recentred: a genuine kernel vector sits
// within about 100 times the truncation floor rho^N, while the smallest NON-kernel singular value
// sits 1e4 or more above it.  A cutoff of 1000 * rho^N therefore separates them; the earlier 1e6
// was loose enough to swallow the first non-kernel value and report one dimension too many at
// weight 2.  The absolute floor only guards against asking for more than the working precision.
function ZeroLevel(rho, nterms : prec := 100)
    return Max(10^(-(prec-30)), 1000 * rho^nterms);
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

// Magma's UnitDisc defaults to Center := (9/10)*i, a point of the upper half plane fixed in advance
// and unrelated to the group.  The series is expanded about that point, so the circumradius rho of
// the fundamental domain -- which fixes both the number of terms (rho^N) and the conditioning of
// the linear system -- is an accident of that default.  At D = 10 it gives a domain spanning
// |w| = 0.05 to 0.82; recentring on the domain's own centroid gives 0.49 to 0.65 and rho 0.648.
// Returns a centre in the upper half plane, and the disc and domain computed there.
function RecentredDisc(G, H, PREC, iters)
    ctr := H ! (9/10*H.1);
    DD := UnitDisc( : Center := ctr, Precision := PREC);
    FD := FundamentalDomain(G, DD);
    for it in [1..iters] do
        vs := [ComplexValue(v) : v in FD];
        c := (&+vs) / #vs;
        if Abs(c) lt 10^-12 then break; end if;
        ctr2 := DiscToPlane(H, DD ! c);
        DD2 := UnitDisc( : Center := ctr2, Precision := PREC);
        FD2 := FundamentalDomain(G, DD2);
        r1 := Max([Abs(ComplexValue(v)) : v in FD]);
        r2 := Max([Abs(ComplexValue(v)) : v in FD2]);
        if r2 ge r1 then break; end if;          // keep the better of the two
        ctr := ctr2; DD := DD2; FD := FD2;
    end for;
    return ctr, DD, FD;
end function;

// Dimension from the largest RELATIVE gap, with no absolute cutoff.
//
// ZeroLevel above assumes the kernel sits near the truncation floor rho^N, which holds when the
// rows come from Gamma alone.  It does NOT hold once Atkin-Lehner rows are included: those are
// limited by the quality of the reduction, not by truncation, so the kernel floats far above
// rho^N.  At D = 6, weight 24, the two kernel values are 7.5e-10 and 1.4e-8 against a truncation
// floor of 1e-16, while the first non-kernel value is 0.0195 -- a ratio of 1.4e6, unmistakable,
// but every kernel value is above the absolute cutoff and the count comes out 0.
// So: call a singular value plausibly zero when it is below relcut times the LARGEST one, and take
// the largest ratio jump among those.  Returns the dimension, the gap ratio, and the two values
// bracketing the gap so a caller can report how unambiguous it was.
function GapDimRel(sv, relcut)
    top := sv[#sv];
    nz := #[s : s in sv | s lt relcut*top];
    if nz eq 0 then return 0, 0.0, 0.0, sv[1]; end if;
    best := 0.0; bi := 0;
    for i in [1..Min(nz, #sv-1)] do
        if sv[i] gt 0 and sv[i+1]/sv[i] gt best then best := sv[i+1]/sv[i]; bi := i; end if;
    end for;
    return bi, best, sv[bi], sv[bi+1];
end function;

// GapDimRel always names a dimension: it returns the index of the largest relative gap even when
// the spectrum is a smooth continuum with no gap at all, which is what happens at D = 15.  This
// one refuses to answer in that case, returning -1.  A dimension is reported only when the
// gap is at least mingap and the first value above it is a real singular value rather than more
// noise, taken as at least 1e-6 of the largest.  The kernel values are NOT required to sit at the
// truncation floor: once the rows come from a reduction rather than from truncation alone they
// float above it, by a factor of 37 at D = 15.  Voight-Willis call a gap of 4
// unsatisfactory in their Example 1 and a gap above 1e19 convincing in their Example 4.
function GapDimStrict(sv, floor, mingap)
    if #sv lt 2 or sv[#sv] le 0 then return -1, 0.0; end if;
    best := 0.0; bi := 0;
    for i in [1..#sv-1] do
        if sv[i] gt 0 and sv[i+1]/sv[i] gt best then best := sv[i+1]/sv[i]; bi := i; end if;
    end for;
    if bi eq 0 or best lt mingap then return -1, best; end if;
    if sv[bi+1] lt 10^-6 * sv[#sv] then return -1, best; end if;
    return bi, best;
end function;
