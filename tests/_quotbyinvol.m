// tests/_quotbyinvol.m -- derive the quotient of a hyperelliptic curve by an Atkin-Lehner
// involution.  Leading underscore: a shared helper, NOT a test (see run_tests.m).
// Its consumer and its negative controls are tests/QuotientByInvolution.m.
//
// WHY THIS EXISTS.  Published sources state an equation for the TOP curve X_0(D,N) and the action
// of each w_m on it, and nothing else; every quotient key is then unpublished.  That is why 10 of
// the 50 bases with an X0_ test check 1 of their 15 committed cover keys.  Quotienting the
// PUBLISHED top by the PUBLISHED involutions produces expected curves that are still external
// data, so this helper EXTENDS an oracle down the subgroup lattice rather than re-deriving
// pipeline output from pipeline output.
//
// ⚠⚠ NEVER FEED IT A COMMITTED MODEL AND CALL THE RESULT AN ORACLE.  The inputs must be
// cover_data[{1}] and ws_data[{1}] -- both external.  Quotienting a committed model by the
// pipeline's own `ws` is circular and proves nothing.
//
// ⚠ IT IS NOT A MODEL-BUILDING LEVER.  Quotienting goes DOWN the subgroup lattice and the model
// gaps are UP it: measured 2026-09-24, this reaches 16 of the 353 unbuilt keys (5%).  The frontier
// is the top curve, which only the Borcherds pipeline builds.
//
// THE METHOD, and why it needs no case analysis.  For sigma a Mobius involution on x:
//   * the invariant is  u = x + sigma(x), falling back to  u = x*sigma(x)  when the sum is
//     constant (the x -> -x case).  Both cannot degenerate: they are the coefficients of
//     (T-x)(T-sigma x), and if both were constant x would be algebraic over Q.
//   * s = x - sigma(x) is ALWAYS anti-invariant, which is why "sum or product" is the only split.
//   * the y-branch is not a case split either.  Write w(y) = phi*y with phi = e/(c*x+d)^(g+1);
//     phi has norm 1 down to Q(u), so by Hilbert 90 there is t with t/sigma(t) = phi, and
//     v = t*y is invariant.  t = r + phi*sigma(r) for almost any r.  phi = 1 gives t = 1, i.e.
//     "y descends"; otherwise t is the anti-invariant y has to be paired with.
//   * v^2 = t^2 f(x) is then sigma-invariant, hence a rational function of u, recovered by linear
//     algebra and CERTIFIED by substitution before it is returned.

// --------------------------------------------------------------------------------------------
// Atkin-Lehner composition.  ⚠ w_m * w_n is w_{m n / gcd(m,n)^2}, NOT w_{mn}.  Multiplying the
// labels instead cost a wrong table on 2026-09-24 (M2*M3*M39 at 39_2 is w_26, not w_78).
function ALCompose(m, n)
    g := GCD(m, n);
    return (m * n) div (g * g);
end function;

// Two matrices describe the SAME map on P(1,w,1) iff they differ by the weighted scaling
// (x,y,z) -> (lam*x, lam^w*y, lam*z).  Used to check a generated group for consistency.
function ALSameMap(M, Mp, w)
    for i, j in [1, 3] do
        if M[i,j] eq 0 and Mp[i,j] ne 0 then return false; end if;
        if M[i,j] ne 0 and Mp[i,j] eq 0 then return false; end if;
    end for;
    lam := 0;
    for i, j in [1, 3] do
        if M[i,j] ne 0 then lam := Mp[i,j] / M[i,j]; break; end if;
    end for;
    if lam eq 0 then return false; end if;
    for i, j in [1, 3] do
        if Mp[i,j] ne lam * M[i,j] then return false; end if;
    end for;
    return Mp[2,2] eq lam^w * M[2,2];
end function;

// Close a set of labelled generator matrices into the full Atkin-Lehner group.
// `gens` is an associative array  m -> 3x3 matrix.  Returns  m -> matrix  for EVERY element.
// ⚠ The labels of the products are computed from the labels of the generators, so a mislabelled
// GENERATOR propagates.  What this does catch is an INCONSISTENT set: if two different words spell
// the same label with genuinely different maps, it says so.
function ALGroupFromGenerators(gens, w)
    all := AssociativeArray();
    ok := true;  why := "";
    for m in Keys(gens) do all[m] := gens[m]; end for;
    repeat
        added := false;
        for m in Keys(all) do
            for n in Keys(gens) do
                p := ALCompose(m, n);
                if p eq 1 then continue; end if;            // w_m * w_m = identity
                P := all[m] * gens[n];
                if IsDefined(all, p) then
                    if not ALSameMap(all[p], P, w) then
                        ok := false;
                        why := Sprintf("two words for w_%o give different maps", p);
                    end if;
                else
                    all[p] := P;  added := true;
                end if;
            end for;
        end for;
    until not added;
    return all, ok, why;
end function;

// --------------------------------------------------------------------------------------------
// Weight-respecting linear maps on P(1,g+1,1) send x,z among themselves and scale y, so with
// X = x/z:   sigma(X) = (a*X+b)/(c*X+d)   and   w(y_affine) = e*y_affine/(c*X+d)^(g+1).
// Row-vector convention throughout: (x,y,z) |-> (x,y,z)*M, matching every ws_data in tests/.
function SigmaData(M)
    return M[1,1], M[3,1], M[1,3], M[3,3], M[2,2];
end function;

// Express a sigma-invariant F in Q(X) as A(u)/B(u).  Linear in the unknown coefficients once
// denominators are cleared; the candidate is certified by substitution before being returned.
function ExpressInU(F, u)
    Nu := Numerator(u);   Du := Denominator(u);
    NF := Numerator(F);   DF := Denominator(F);
    n := Max([Degree(NF), Degree(DF)]);
    for k in [Maximum(1, n div 2) .. n + 2] do          // grow the ansatz until it solves
        cols := [];
        for j in [0..k] do Append(~cols,  NF * Nu^j * Du^(k-j)); end for;    // beta_j
        for i in [0..k] do Append(~cols, -DF * Nu^i * Du^(k-i)); end for;    // alpha_i
        deg := Max([Degree(c) : c in cols]);
        Mat := Matrix(Rationals(), #cols, deg+1,
                      [[Coefficient(c, r) : r in [0..deg]] : c in cols]);
        for v in Basis(Kernel(Mat)) do
            s := Eltseq(v);
            T := PolynomialRing(Rationals());
            B := T ! s[1..k+1];
            A := T ! s[k+2..2*k+2];
            if B eq 0 then continue; end if;
            gg := GCD(A, B);
            if Degree(gg) gt 0 then A := A div gg; B := B div gg; end if;
            if B ne 0 and Evaluate(A, u)/Evaluate(B, u) eq F then return true, A, B; end if;
        end for;
    end for;
    return false, 0, 0;
end function;

// --------------------------------------------------------------------------------------------
// The quotient of  y^2 = f(x)  by the involution with matrix M.
// Returns  ok, Fq, note.  Fq = 0 with ok means the quotient is P^1 (the x-line).
function QuotientByInvolution(f, M)
    Pol<T> := PolynomialRing(Rationals());
    f := Pol ! f;
    g  := (Degree(f) - 1) div 2;                 // deg f = 2g+1 or 2g+2
    wy := g + 1;                                 // weight of y

    a, b, c, d, e := SigmaData(M);

    FF<X> := FunctionField(Rationals());
    den := c*X + d;
    if den eq 0 then return false, 0, "degenerate: c*X+d vanishes"; end if;
    sig := (a*X + b) / den;

    // CONTROL: sigma must be an involution.
    if Evaluate(sig, sig) ne X then return false, 0, "sigma is not an involution"; end if;

    // sigma = id is the hyperelliptic involution (x,y) -> (x,-y); the quotient is the x-line.
    if sig eq X then
        if e eq -1 then return true, Pol!0, "P1 (hyperelliptic involution)"; end if;
        return false, 0, "sigma is the identity and e /= -1: this is the identity map";
    end if;

    // CONTROL, the load-bearing one: the matrix must PRESERVE the curve.  y_new^2 = f(X_new)
    // means e^2 f(X) = f(sigma X) * (c X + d)^(2g+2).  Without this a wrong matrix still yields a
    // plausible-looking curve -- the exact failure mode CLAUDE.md warns about.
    if e^2 * Evaluate(FF!f, X) ne Evaluate(FF!f, sig) * den^(2*g+2) then
        return false, 0, "matrix does NOT preserve the curve";
    end if;

    u := X + sig;  note := "u = x + sigma(x)";
    if u in Rationals() then
        u := X * sig;  note := "u = x * sigma(x) (sum degenerate)";
        if u in Rationals() then return false, 0, "both x+sigma and x*sigma are constant"; end if;
    end if;

    phi := e / den^wy;
    if phi * Evaluate(phi, sig) ne 1 then return false, 0, "phi has norm /= 1"; end if;
    t := FF ! 0;
    if phi eq 1 then
        t := FF ! 1;  note cat:= " ; y descends";
    else
        for r in [FF!1, X, X^2, X+1, X^2+X+1] do
            cand := r + phi * Evaluate(r, sig);
            if cand ne 0 then t := cand; break; end if;
        end for;
        if t eq 0 then return false, 0, "Hilbert 90: no t found"; end if;
        note cat:= " ; y paired with an anti-invariant";
    end if;
    if Evaluate(t, sig) * phi ne t then return false, 0, "Hilbert 90 certificate failed"; end if;

    Fv := t^2 * Evaluate(FF!f, X);
    if Evaluate(Fv, sig) ne Fv then return false, 0, "t^2 f is not sigma-invariant"; end if;
    ok, A, B := ExpressInU(Fv, u);
    if not ok then return false, 0, "could not express t^2 f as a rational function of u"; end if;

    // v^2 = A(u)/B(u)  =>  (B v)^2 = A B.  y^2 = F and y^2 = F/s^2 are the same curve, so drop
    // every squared factor but KEEP the constant: its square class is a real invariant.
    Fq := A * B;
    if Fq eq 0 then return false, 0, "quotient form vanishes"; end if;
    red := Pol ! LeadingCoefficient(Fq);
    for pr in Factorisation(Fq) do
        if IsOdd(pr[2]) then red *:= pr[1]; end if;
    end for;
    return true, red, note;
end function;

// --------------------------------------------------------------------------------------------
// ⚠ THE DERIVED MODEL IS ISOMORPHIC TO A COMMITTED ONE, NOT EQUAL TO IT.  At 14_1, w_2 gives
// -u^2-13u-128 against a committed -u^2+13u-128: u = x*sigma(x) = -x^2 differs from the committed
// normalisation by u -> -u.  Compare by isomorphism class, never coefficient-wise.
// Genus 0 must NOT go through IsIsomorphic (Magma 2.29-10 regression, see BorcherdsProducts.m);
// compare the Brauer class exactly as tests/ConicClasses.m does.
function QuotientConicClass(f)
    if Degree(f) ne 2 then return false, []; end if;
    a := Coefficient(f, 2);
    disc := Coefficient(f, 1)^2 - 4*a*Coefficient(f, 0);
    if a eq 0 or disc eq 0 then return false, []; end if;
    return true, Sort(RamifiedPrimes(QuaternionAlgebra<Rationals() | a, disc>));
end function;

// Compare a derived quotient against a committed/expected curve.  Returns  same, why.
function QuotientMatches(Fder, Cexp)
    fe := HyperellipticPolynomials(SimplifiedModel(Cexp));
    gE := Genus(Cexp);
    if Fder eq 0 then
        return gE eq 0 and Degree(fe) le 2, Sprintf("P1 vs genus %o", gE);
    end if;
    Cder := HyperellipticCurve(Fder);
    gD := Genus(Cder);
    if gD ne gE then return false, Sprintf("genus %o vs %o", gD, gE); end if;
    if gD ge 1 then
        iso := IsIsomorphic(Cder, Cexp);
        return iso, Sprintf("genus %o, IsIsomorphic=%o", gD, iso);
    end if;
    okd, rd := QuotientConicClass(Fder);
    oke, re := QuotientConicClass(fe);
    if not okd or not oke then
        return Degree(Fder) le 2 and Degree(fe) le 2, "both split genus 0";
    end if;
    return rd eq re, Sprintf("conic class %o vs %o", rd, re);
end function;
