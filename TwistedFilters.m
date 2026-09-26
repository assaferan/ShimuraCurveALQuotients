// Twisted-trace and twisted-Weil-polynomial non-hyperellipticity filters.
//
// C = X_0^D(N)/W, h an involution of C defined over Q, p a prime not dividing DN.  The quadratic
// twist C_h of C by h has Frob_p acting on H^1(C_h) as Frob_p o h on H^1(C).  If C is hyperelliptic
// then so is C_h (the hyperelliptic involution is central in Aut(C)), which gives two tests:
//
//  * FilterByTwistedTrace.  #C_h(F_q) <= 2(q+1), i.e. for q = p^v,
//        tr := Tr((T_{p^v} - p T_{p^{v-2}}) o h | S_2(DN; W = sigma)^{D-new}) >= -(q+1).
//    |tr| <= 2g sqrt(q), so a violation needs q < 4g^2, the range FilterByTrace uses.
//  * FilterByTwistedWeilPolynomial.  P_{C_h}(t) = P_+(t) P_-(-t), with P_+- the product of
//    (t^2 - a t + p) over the T_p-eigenvalues a on the h = +-1 eigenspaces, must be the Weil
//    polynomial of a hyperelliptic curve over F_p: a line of data/hypg<g>q<p>.txt.  Those LMFDB
//    tables are complete for g = 3, p <= 23; g = 4, p <= 5; g = 5, 6, p = 2, and only those are used.
//    Parity and 2-rank criteria add nothing, since P_{C_h} = P_C (mod 2).
//
// Here w_m acts by sigma_m w_m, sigma_m = (-1)^omega(gcd(m, D)), on the D-new cuspidal modular
// symbols of level DN (sign 0, so each form appears twice), and C's H^1 is the W-fixed part K.
// The involutions h used are:
//   * the residual Atkin-Lehner w_Q, Q notin W;
//   * S2, V2, V3, their pairwise products, and these times any w_Q, under the descent conditions of
//     CheckModularNonALInvolutionModSym (S2 needs every w in W odd, and only odd Q);
//   * V3 only when 9 in W: Galois acts by sigma(V3) = V3 W9, so V3 is defined over Q on C only then.
// An op is kept only if it preserves K, is an involution on K (M^2 = 1), and commutes with T_p on K
// for every prime p the test uses on that curve (a safety net: an op defined over Q commutes with
// every T_p).  h = 1 is left out: it is FilterByTrace / FilterByWeilPolynomial.
//
// Why the results are trustworthy (prototypes: sweeps/twisted_trace/twist.m and
// sweeps/twisted_weil/tw.m on branch twisted-trace-sweep, with README): the sweeps found 0
// violations on every recorded-hyperelliptic control curve; the AL-twisted traces agree with
// Eichler-Selberg (TraceDNewALFixed); the reviewer's independent code reproduced the trace hits; and
// the Weil code at h = 1 reproduces all 60 recorded "WeilPolynomial with p = x" verdicts.
// tests/TwistedTraceFilter.m and tests/TwistedWeilFilter.m pin the known values.

intrinsic TwistedModSymMaxLevel() -> RngIntElt
{Maximum level D*N at which FilterByTwistedTrace and FilterByTwistedWeilPolynomial run; curves at
larger levels are left undetermined.  Both need the modular symbols space of level D*N and its Hecke
operators, which at D*N of several thousand costs hours per level.  The default is above every level
in the current curve list (at most 15330), so nothing is skipped; lower it to bound the cost.}
    return 20000;
end intrinsic;

// Primes whose Weil polynomial tables are complete, by genus.
function twTablePrimes(g)
    if g eq 3 then return [2,3,5,7,11,13,17,19,23]; end if;
    if g eq 4 then return [2,3,5]; end if;
    if g in {5,6} then return [2]; end if;
    return [];
end function;

// Primes used by the twisted trace on a curve of genus g at level L: p good, p < 4g^2.
function twTracePrimes(g, L)
    return [p : p in PrimesUpTo(4*g^2 - 1) | L mod p ne 0];
end function;

function twWeilPrimes(g, L)
    return [p : p in twTablePrimes(g) | L mod p ne 0];
end function;

// Level data at (D,N): the D-new cuspidal modular symbols (basis B), sigma-twisted AL matrices,
// the non-AL S2, V2, V3 that exist at this level, and T_p for p in ps, all on that basis.
function twLevelData(D, N, ps)
    L := D*N;
    MDN := ModularSymbols(L, 2, 0);
    SDN := CuspidalSubspace(MDN);
    for p in PrimeDivisors(D) do SDN := NewSubspace(SDN, p); end for;
    B := Matrix([Representation(v) : v in Basis(SDN)]);
    n := Nrows(B);
    MA := MatrixAlgebra(Rationals(), n);
    ALm := AssociativeArray();
    for Q in [Q : Q in Divisors(L) | GCD(Q, L div Q) eq 1] do
        chi := (-1)^#PrimeDivisors(GCD(Q, D));
        ALm[Q] := MA!(chi*Solution(B, B*AtkinLehnerOperator(MDN, Q)));
    end for;
    Vfull := AssociativeArray();
    if N mod 4 eq 0 then Vfull["S2"] := MA!ModularNonALOperatorOnSubspace(2, N, B, MDN, true); end if;
    if N mod 8 eq 0 then Vfull["V2"] := MA!ModularNonALOperatorOnSubspace(2, N, B, MDN, false); end if;
    if Valuation(N, 3) eq 2 then Vfull["V3"] := MA!ModularNonALOperatorOnSubspace(3, N, B, MDN, false); end if;
    Tp := AssociativeArray();
    for p in ps do Tp[p] := MA!Solution(B, B*HeckeOperator(MDN, p)); end for;
    vprintf ShimuraQuotients, 2: "twisted: level (%o,%o), D-new dimension %o, %o Hecke operators\n", D, N, n, #ps;
    return <D, N, n, ALm, Vfull, Tp>;
end function;

// The data of curve X on the level data LD: the restrictions to K = H^1(X) of T_p (p in ps) and of
// the admissible involutions h != 1, as a list of <name, matrix>.
function twCurveData(X, LD, ps)
    D, N, n, ALm, Vfull, Tp := Explode(LD);
    W := X`W; L := D*N;
    K := VectorSpace(Rationals(), n);
    for w in W do
        if w ne 1 then K meet:= Kernel(ALm[w] - 1); end if;
    end for;
    BK := BasisMatrix(K);
    dK := Nrows(BK);
    error if dK ne 2*X`g, Sprintf("twisted: W-fixed D-new dimension %o is not 2g = %o on %o", dK, 2*X`g, X);
    MK := MatrixAlgebra(Rationals(), dK);
    IK := MK!1;
    Tk := AssociativeArray();
    for p in ps do
        ok, S := IsConsistent(BK, BK*Tp[p]);
        assert ok;       // T_p commutes with the ALs, so it preserves K
        Tk[p] := MK!S;
    end for;
    als := Sort([Q : Q in Keys(ALm)]);
    cand := [<"w" cat IntegerToString(Q), ALm[Q]> : Q in als | Q notin W];
    vn := [];
    if IsDefined(Vfull, "S2") and &and[IsOdd(w) : w in W] then Append(~vn, "S2"); end if;
    if IsDefined(Vfull, "V2") then Append(~vn, "V2"); end if;
    if IsDefined(Vfull, "V3") and (9 in W) then Append(~vn, "V3"); end if;   // Q-rational V3 only
    allv := [<a, Vfull[a]> : a in vn] cat [<a cat "*" cat b, Vfull[a]*Vfull[b]> : a, b in vn | a ne b];
    for vv in allv do
        for Q in als do
            if (Q in W) and (Q ne 1) then continue; end if;
            if ("S2" in vv[1]) and IsEven(Q) then continue; end if;
            Append(~cand, <Q eq 1 select vv[1] else vv[1] cat "*w" cat IntegerToString(Q), vv[2]*ALm[Q]>);
        end for;
    end for;
    ops := [];
    for o in cand do
        ok, S := IsConsistent(BK, BK*o[2]);
        if not ok then continue; end if;            // does not descend to K
        M := MK!S;
        if M eq IK then continue; end if;           // h = 1 on C
        if M^2 ne IK then
            vprintf ShimuraQuotients, 2: "twisted: %o is not an involution on %o\n", o[1], X;
            continue;
        end if;
        if exists{p : p in ps | Tk[p]*M ne M*Tk[p]} then
            vprintf ShimuraQuotients, 1: "twisted: %o does not commute with every T_p on %o; dropped\n", o[1], X;
            continue;
        end if;
        if exists{x : x in ops | x[2] eq M} then continue; end if;
        Append(~ops, <o[1], M>);
    end for;
    return Tk, ops, IK;
end function;

// First twisted-trace violation on X, over p in ps ascending, then v, then the ops.
function twTraceCheck(X, LD)
    g := X`g;
    ps := twTracePrimes(g, X`D*X`N);
    Tk, ops, IK := twCurveData(X, LD, ps);
    QM := 4*g^2 - 1;
    for p in ps do
        T := Tk[p];
        s0 := 2*IK; s1 := T; v := 1;       // s_v = T_{p^v} - p T_{p^{v-2}} = alpha^v + beta^v
        while p^v le QM do
            q := p^v;
            for o in ops do
                tr := Integers()!(Trace(s1*o[2]) / 2);
                // #C_h(F_q) = q + 1 - tr >= 0 on any curve
                error if tr gt q + 1, Sprintf("twisted trace %o > q + 1 = %o on %o, h = %o", tr, q + 1, X, o[1]);
                if tr lt -(q + 1) then return false, o[1], p, v, tr; end if;
            end for;
            s2 := T*s1 - p*s0; s0 := s1; s1 := s2; v +:= 1;
        end while;
    end for;
    return true, _, _, _, _;
end function;

// prod (t^2 - a t + p) over the roots a of c, where chiT = c^2 is the charpoly of T_p on a
// (sign 0) modular symbols subspace: t^deg(c) c(t + p/t).
function twWeilFactor(chiT, p)
    R2 := PolynomialRing(Rationals());
    fa := Factorization(R2!chiT);
    assert &and[e[2] mod 2 eq 0 : e in fa];
    c := &*[R2 | e[1]^(e[2] div 2) : e in fa];
    d := Degree(c); cf := Coefficients(c);
    Zt<u> := PolynomialRing(Integers());
    return &+[Zt | Integers()!cf[i+1] * (u^2 + p)^i * u^(d-i) : i in [0..d]];
end function;

// P_h(t) = P_+(t) P_-(-t) for the involution M (or the identity) on K, with sanity checks.
function twWeilPolynomial(T, M, IK, p, g)
    Zt<u> := PolynomialRing(Integers());
    Kp := Kernel(M - IK); Km := Kernel(M + IK);
    Pp := Zt!1; Pm := Zt!1;
    if Dimension(Kp) gt 0 then
        Bp := BasisMatrix(Kp); Pp := twWeilFactor(CharacteristicPolynomial(Solution(Bp, Bp*T)), p);
    end if;
    if Dimension(Km) gt 0 then
        Bm := BasisMatrix(Km); Pm := twWeilFactor(CharacteristicPolynomial(Solution(Bm, Bm*T)), p);
    end if;
    Ph := Pp * Evaluate(Pm, -u);
    assert Degree(Ph) eq 2*g and LeadingCoefficient(Ph) eq 1 and Coefficient(Ph, 0) eq p^g;
    assert Coefficient(Ph, 2*g-1) eq -Trace(T*M)/2;
    return Ph;
end function;

// The hyperelliptic Weil polynomials of genus g over F_p, as the strings "[1,a1,...,p^g]" of
// data/hypg<g>q<p>.txt (the lines LMFDBweilpolys evaluates; a set of strings loads far faster).
function twWeilTable(g, p)
    return {l : l in Split(Read(Sprintf("data/hypg%oq%o.txt", g, p)), "\n") | #l gt 0};
end function;

function twWeilKey(P)
    return "[" cat Join([IntegerToString(x) : x in Reverse(Coefficients(P))], ",") cat "]";
end function;

// First twisted Weil polynomial of X not in the tables, over p in ps ascending, then the ops.
function twWeilCheck(X, LD, TAB)
    g := X`g;
    ps := twWeilPrimes(g, X`D*X`N);
    if #ps eq 0 then return true, _, _; end if;
    Tk, ops, IK := twCurveData(X, LD, ps);
    for p in ps do
        P1 := twWeilPolynomial(Tk[p], IK, IK, p, g);
        for o in ops do
            Ph := twWeilPolynomial(Tk[p], o[2], IK, p, g);
            assert ChangeRing(Ph, GF(2)) eq ChangeRing(P1, GF(2));
            if twWeilKey(Ph) notin TAB[<g, p>] then return false, o[1], p; end if;
        end for;
    end for;
    return true, _, _;
end function;

function twLoadTables(gs)
    TAB := AssociativeArray();
    for g in gs do
        for p in twTablePrimes(g) do TAB[<g, p>] := twWeilTable(g, p); end for;
    end for;
    return TAB;
end function;

intrinsic CheckTwistedTrace(X::ShimuraQuot) -> BoolElt, MonStgElt, RngIntElt, RngIntElt, RngIntElt
{Returns false if some admissible involution h != 1 of X defined over Q and some q = p^v < 4g^2,
p good, have Tr((T_(p^v) - p T_(p^(v-2))) o h) < -(q+1), proving X non-hyperelliptic; then also
returns the name of h, p, v and the trace.  Returns true otherwise.}
    assert X`g ge 2;
    LD := twLevelData(X`D, X`N, twTracePrimes(X`g, X`D*X`N));
    ok, h, p, v, tr := twTraceCheck(X, LD);
    if ok then return true, _, _, _, _; end if;
    return false, h, p, v, tr;
end intrinsic;

intrinsic CheckTwistedWeilPolynomial(X::ShimuraQuot) -> BoolElt, MonStgElt, RngIntElt
{Returns false if for some admissible involution h != 1 of X defined over Q and some table prime p
the twisted Weil polynomial P_(X_h) is not the Weil polynomial of a hyperelliptic curve over F_p,
proving X non-hyperelliptic; then also returns the name of h and p.  Returns true otherwise.}
    if #twWeilPrimes(X`g, X`D*X`N) eq 0 then return true, _, _; end if;
    LD := twLevelData(X`D, X`N, twWeilPrimes(X`g, X`D*X`N));
    ok, h, p := twWeilCheck(X, LD, twLoadTables({X`g}));
    if ok then return true, _, _; end if;
    return false, h, p;
end intrinsic;

intrinsic TwistedWeilPolynomials(X::ShimuraQuot, p::RngIntElt) -> SeqEnum
{The twisted Weil polynomials P_(X_h) at the good prime p, as a list of <name of h, P_h>, for h = 1
followed by the admissible involutions h != 1 of X defined over Q.  At h = 1 this is the Weil
polynomial of X.}
    require (X`D*X`N) mod p ne 0 : "p must be a good prime";
    LD := twLevelData(X`D, X`N, [p]);
    Tk, ops, IK := twCurveData(X, LD, [p]);
    return [<"1", twWeilPolynomial(Tk[p], IK, IK, p, X`g)>] cat
           [<o[1], twWeilPolynomial(Tk[p], o[2], IK, p, X`g)> : o in ops];
end intrinsic;

// The undecided curves of genus >= 3 in curves, grouped by level (D,N) in increasing D*N; only the
// curves with some prime to test (per primes_of(g, D*N)) and at levels within the cap.
function twLevels(curves, primes_of)
    levels := AssociativeArray();
    for i->X in curves do
        if assigned X`IsSubhyp then continue; end if;
        if X`g lt 3 then continue; end if;
        if #primes_of(X`g, X`D*X`N) eq 0 then continue; end if;
        if X`D*X`N gt TwistedModSymMaxLevel() then
            vprintf ShimuraQuotients, 2: "twisted: %o exceeds the level cap\n", X;
            continue;
        end if;
        key := <X`D, X`N>;
        if not IsDefined(levels, key) then levels[key] := []; end if;
        Append(~levels[key], i);
    end for;
    keys := [k : k in Keys(levels)];
    Sort(~keys, func<a, b | a[1]*a[2] - b[1]*b[2]>);
    return keys, levels;
end function;

intrinsic FilterByTwistedTrace(~curves::SeqEnum)
{Mark as non-hyperelliptic the curves that fail the twisted trace test (see CheckTwistedTrace).
The modular symbols are computed once per level.}
    keys, levels := twLevels(curves, twTracePrimes);
    for j->key in keys do
        vprintf ShimuraQuotients, 1: "twisted trace: level %o/%o, (D,N) = %o, %o curves\n", j, #keys, key, #levels[key];
        idx := levels[key];
        ps := twTracePrimes(Max([curves[i]`g : i in idx]), key[1]*key[2]);
        LD := twLevelData(key[1], key[2], ps);
        for i in idx do
            ok, h, p, v, tr := twTraceCheck(curves[i], LD);
            if not ok then
                vprintf ShimuraQuotients, 2: "curve %o: h = %o, q = %o^%o, trace %o\n", curves[i]`CurveID, h, p, v, tr;
                curves[i]`IsSubhyp := false;
                curves[i]`IsHyp := false;
                curves[i]`TestInWhichProved := Sprintf("TwistedTrace, h = %o with p^v = %o^%o", h, p, v);
            end if;
        end for;
    end for;
end intrinsic;

intrinsic FilterByTwistedWeilPolynomial(~curves::SeqEnum)
{Mark as non-hyperelliptic the curves whose twisted Weil polynomials fail the LMFDB hyperelliptic
tables (see CheckTwistedWeilPolynomial).  The modular symbols are computed once per level.}
    keys, levels := twLevels(curves, twWeilPrimes);
    if #keys eq 0 then return; end if;
    TAB := twLoadTables(&join[{curves[i]`g : i in levels[key]} : key in keys]);
    for j->key in keys do
        vprintf ShimuraQuotients, 1: "twisted Weil: level %o/%o, (D,N) = %o, %o curves\n", j, #keys, key, #levels[key];
        idx := levels[key];
        ps := &join[Seqset(twWeilPrimes(curves[i]`g, key[1]*key[2])) : i in idx];
        LD := twLevelData(key[1], key[2], Sort(SetToSequence(ps)));
        for i in idx do
            ok, h, p := twWeilCheck(curves[i], LD, TAB);
            if not ok then
                vprintf ShimuraQuotients, 2: "curve %o: h = %o, p = %o\n", curves[i]`CurveID, h, p;
                curves[i]`IsSubhyp := false;
                curves[i]`IsHyp := false;
                curves[i]`TestInWhichProved := Sprintf("TwistedWeilPolynomial h = %o with p = %o", h, p);
            end if;
        end for;
    end for;
end intrinsic;
