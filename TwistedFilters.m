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

import !"Geometry/ModSym/operators.m" : ActionOnModularSymbolsBasis, Heilbronn, TnSparse;

// There is deliberately no level cap: every level is run, and the largest (D*N up to 15330) take
// hours each.  See docs/RUNNING_PIPELINE.md.

// Can X carry an admissible involution h != 1 at all?  (A residual AL, S2 with every w in W odd,
// V2, or V3 with 9 in W.)  If not, the twisted tests have nothing to do and the modular symbols of
// the level are not built.  On a star curve (W the full AL group) only V2 and V3 can qualify.
function twHasOps(X)
    N := X`N; W := X`W; L := X`D*N;
    if exists{Q : Q in Divisors(L) | GCD(Q, L div Q) eq 1 and Q notin W} then return true; end if;
    if (N mod 4 eq 0) and &and[IsOdd(w) : w in W] then return true; end if;
    if N mod 8 eq 0 then return true; end if;
    return (Valuation(N, 3) eq 2) and (9 in W);
end function;

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

// Pivot columns of a matrix in echelon form.
function twPivots(E)
    return [Min([j : j in [1..Ncols(E)] | E[i,j] ne 0]) : i in [1..Nrows(E)]];
end function;

// Level data at (D,N).  Everything is kept at the size of the D-new space or of K = H^1(X), never
// as a full operator on the D-new space built by Solution against its (dense, large-height) basis,
// which is what made the large levels slow:
//  * B and Phi, the echelonized bases of the D-new cuspidal modular symbols SDN and of its dual
//    (the functionals vanishing on the Hecke-stable complement: Eisenstein and D-old).  For an
//    echelon basis the coordinates of a vector are its entries at the pivot columns, so an operator
//    A preserving SDN is B * A[:, pivots] on it, with no Solution.
//  * the sigma-twisted AL matrices only for the prime powers l^a || D*N, on B and on Phi; w_Q for
//    composite Q is their product (W_a W_b = W_ab exactly on modular symbols of weight 2, trivial
//    character).  Before this every Hall divisor Q got its own matrix (64 at level 8190).
//  * S_mu, S_mu^-1 and W_(mu^v) on the ambient, for the V_mu = S_mu W S_mu^-1 of the level.
//    S_mu^-1 is the action of [mu,-1,0,mu] (= mu^2 S_mu^-1, and scalars act trivially in weight 2),
//    so nothing is inverted; these are applied to the few rows of K only (twCurveData).
//  * T_p is not computed here at all: twCurveData gets it on K from d = dim K Manin symbols.
function twLevelData(D, N)
    L := D*N;
    MDN := ModularSymbols(L, 2, 0);
    SDN := CuspidalSubspace(MDN);
    for p in PrimeDivisors(D) do SDN := NewSubspace(SDN, p); end for;
    B := EchelonForm(Matrix([Representation(v) : v in Basis(SDN)]));
    Phi := EchelonForm(BasisMatrix(DualVectorSpace(SDN)));
    pivB := twPivots(B); pivPhi := twPivots(Phi);
    n := Nrows(B);
    ALB := AssociativeArray(); ALPhi := AssociativeArray();
    for l in PrimeDivisors(L) do
        q := l^Valuation(L, l);
        A := ((D mod l eq 0) select -1 else 1) * AtkinLehnerOperator(MDN, q);
        ALB[q] := B * ColumnSubmatrix(A, pivB);                              // A on B
        ALPhi[q] := Phi * Transpose(Matrix([A[i] : i in pivPhi]));           // A^T on Phi
    end for;
    Vamb := AssociativeArray();
    // Only where twCurveData can use them: S2 needs 4 | N (and V2 8 | N), V3 needs v_3(N) = 2.
    for mu in [2, 3] do
        if (mu eq 2) and (N mod 4 ne 0) then continue; end if;
        if (mu eq 3) and (Valuation(N, 3) ne 2) then continue; end if;
        Vamb[mu] := <ActionOnModularSymbolsBasis([mu,1,0,mu], MDN), AtkinLehnerOperator(MDN, mu^Valuation(N, mu)),
                     ActionOnModularSymbolsBasis([mu,-1,0,mu], MDN)>;
    end for;
    vprintf ShimuraQuotients, 2: "twisted: level (%o,%o), D-new dimension %o, ambient dimension %o\n", D, N, n, Dimension(MDN);
    return <D, N, n, MDN, B, Phi, ALB, ALPhi, Vamb>;
end function;

// The sigma-twisted w_Q (Q a Hall divisor) on a space whose prime-power AL matrices are AL, as the
// product of the prime-power ones.
function twALProduct(AL, Q, L, MA)
    M := MA!1;
    for l in PrimeDivisors(Q) do M := M * AL[l^Valuation(L, l)]; end for;
    return M;
end function;

// Echelon row basis (in the coordinates of the space) of the W-fixed part, cutting by a basis gens
// of W over F_2.  Each w_Q preserves the part fixed by the previous ones (the ALs commute), so its
// matrix there is read off at the pivots.  w_Q is applied to the rows one prime power at a time.
function twFixed(AL, gens, L, n)
    K := MatrixAlgebra(Rationals(), n)!1;
    for Q in gens do
        R := K;
        for l in PrimeDivisors(Q) do R := R * AL[l^Valuation(L, l)]; end for;
        C := ColumnSubmatrix(R, twPivots(K));
        K := EchelonForm(KernelMatrix(C - 1) * K);
    end for;
    return K;
end function;

// Rows R (ambient coordinates) under S2, or V_mu = S_mu W S_mu^-1 (on rows: ((x S) W) S^-1), on the
// ambient.
function twActV(R, a, Vamb)
    if a eq "S2" then return R * Vamb[2][1]; end if;
    t := Vamb[a eq "V2" select 2 else 3];
    return ((R * t[1]) * t[2]) * t[3];
end function;

// The data of curve X on the level data LD: the restrictions to K = H^1(X) of T_p (p in ps) and of
// the admissible involutions h != 1, as a list of <name, matrix>, all in one basis of K.
function twCurveData(X, LD, ps)
    D, N, n, MDN, B, Phi, ALB, ALPhi, Vamb := Explode(LD);
    W := X`W; L := D*N;
    // W is cut out by a basis over F_2 (fixed by it = fixed by all of W, as w_Q, sigma are
    // multiplicative).  Not the prime-power ALs: W need not be generated by them (e.g. <w6>).
    gens := []; span := {1};
    for w in Sort(SetToSequence(W)) do
        if w notin span then Append(~gens, w); span join:= {AtkinLehnerMul(w, s, L) : s in span}; end if;
    end for;
    Kc := twFixed(ALB, gens, L, n);
    BK := Kc * B;                                  // K, in ambient coordinates
    dK := Nrows(BK);
    // Deliberately an error, not a skip: a mismatch means the space is not H^1 of this curve (wrong
    // genus, W or sign convention), and then no verdict of this filter can be trusted.  It stops
    // the (parallel) stage; the prototypes never met it on any curve.
    error if dK ne 2*X`g, Sprintf("twisted filters: BADDIM on curve %o (%o): the W-fixed D-new modular symbols have dimension %o, not 2g = %o; the genus or W is inconsistent, so the stage stops rather than risk a wrong verdict", assigned X`CurveID select X`CurveID else "?", X, dK, 2*X`g);
    PhiK := twFixed(ALPhi, gens, L, n) * Phi;      // functionals vanishing on the complement of K
    assert Nrows(PhiK) eq dK;
    MK := MatrixAlgebra(Rationals(), dK);
    IK := MK!1;
    // T_p on K from dK Manin symbols (Stein's trick, as HeckeOperator does on a subspace, but with
    // PhiK in place of its dense ambient projection matrix): T_p preserves K and its complement, so
    // (x T_p) PhiK^T = (x PhiK^T) [T_p] for every ambient x; take x = the symbols E at the pivots.
    PhiKT := Transpose(PhiK);
    C := MK!(BK * PhiKT);                          // basis change: BK-coordinates -> PhiK-coordinates
    Ci := C^(-1);
    E := twPivots(EchelonForm(PhiK));
    P0i := (MK!Matrix([PhiKT[e] : e in E]))^(-1);
    Tk := AssociativeArray();
    for p in ps do
        H := Heilbronn(MDN, p, false);
        ET := Matrix([TnSparse(MDN, H, [<1, e>]) : e in E]);
        Tk[p] := MK!(C * P0i * (ET * PhiKT) * Ci);
    end for;
    // The same safety net as before (T_p preserves K), now also checking the trick: at the smallest
    // p, where the ambient T_p is cheapest, compare with the ambient operator.  One prime by design:
    // checking every p would bring back the ambient T_p work that the trick exists to avoid.
    if #ps gt 0 then
        ok, S := IsConsistent(BK, BK * HeckeOperator(MDN, ps[1]));
        assert ok and MK!S eq Tk[ps[1]];
    end if;
    // The prime-power ALs on K (they preserve it: the ALs commute), in the BK basis.
    ALK := AssociativeArray();
    pivK := twPivots(Kc);
    for q in Keys(ALB) do ALK[q] := MK!ColumnSubmatrix(Kc * ALB[q], pivK); end for;
    als := [Q : Q in Divisors(L) | GCD(Q, L div Q) eq 1];
    cand := [<"w" cat IntegerToString(Q), twALProduct(ALK, Q, L, MK)> : Q in als | Q notin W];
    vn := [];
    if IsDefined(Vamb, 2) and (N mod 4 eq 0) and &and[IsOdd(w) : w in W] then Append(~vn, "S2"); end if;
    if IsDefined(Vamb, 2) and (N mod 8 eq 0) then Append(~vn, "V2"); end if;
    if IsDefined(Vamb, 3) and (9 in W) then Append(~vn, "V3"); end if;   // Q-rational V3 only
    allv := [<[a], a> : a in vn] cat [<[a, b], a cat "*" cat b> : a, b in vn | a ne b];
    for vv in allv do
        R := BK;
        for a in vv[1] do R := twActV(R, a, Vamb); end for;         // x (V_a V_b) = (x V_a) V_b
        // V w_Q preserves K iff V does (w_Q is invertible on K), so when V does not, none of the
        // V w_Q do either and all of them are skipped below, as before.
        ok, S := IsConsistent(BK, R);
        if not ok then continue; end if;
        MV := MK!S;
        for Q in als do
            if (Q in W) and (Q ne 1) then continue; end if;
            if ("S2" in vv[2]) and IsEven(Q) then continue; end if;
            Append(~cand, <Q eq 1 select vv[2] else vv[2] cat "*w" cat IntegerToString(Q), MV * twALProduct(ALK, Q, L, MK)>);
        end for;
    end for;
    ops := [];
    for o in cand do
        M := o[2];
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

// First twisted-trace violation on X, over the primes ps (each good, < 4g^2) ascending, then v, then
// the ops; Tk, ops, IK are twCurveData's output, on a prime list containing ps.
function twTraceCheckOn(X, Tk, ops, IK, ps)
    g := X`g;
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

// First twisted-trace violation on X, over all the primes the filter uses on X.
function twTraceCheck(X, LD)
    ps := twTracePrimes(X`g, X`D*X`N);
    Tk, ops, IK := twCurveData(X, LD, ps);
    return twTraceCheckOn(X, Tk, ops, IK, ps);
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
    LD := twLevelData(X`D, X`N);
    ok, h, p, v, tr := twTraceCheck(X, LD);
    if ok then return true, _, _, _, _; end if;
    return false, h, p, v, tr;
end intrinsic;

intrinsic CheckTwistedWeilPolynomial(X::ShimuraQuot) -> BoolElt, MonStgElt, RngIntElt
{Returns false if for some admissible involution h != 1 of X defined over Q and some table prime p
the twisted Weil polynomial P_(X_h) is not the Weil polynomial of a hyperelliptic curve over F_p,
proving X non-hyperelliptic; then also returns the name of h and p.  Returns true otherwise.}
    if #twWeilPrimes(X`g, X`D*X`N) eq 0 then return true, _, _; end if;
    LD := twLevelData(X`D, X`N);
    ok, h, p := twWeilCheck(X, LD, twLoadTables({X`g}));
    if ok then return true, _, _; end if;
    return false, h, p;
end intrinsic;

intrinsic CheckTwistedAtPrime(X::ShimuraQuot, p::RngIntElt) -> BoolElt, MonStgElt
{The twisted tests at the single good prime p only.  Returns false and a description if, for some
admissible involution h != 1 of X defined over Q, either the twist X_h has more than 2q+2 points
over F_q for some q = p^v < 4g^2 (the test of CheckTwistedTrace), or p is a table prime and the
Weil polynomial of X_h at p is not that of a hyperelliptic curve over F_p (the test of
CheckTwistedWeilPolynomial).  As X_h is isomorphic to X over the algebraic closure of F_p, X is
then not hyperelliptic over it.  Returns true otherwise; the modular symbols are not built when X
has no admissible h or neither test applies at p.}
    require (X`D*X`N) mod p ne 0 : "p must be a good prime";
    g := X`g;
    assert g ge 2;
    trace_ok := p lt 4*g^2;
    weil_ok := p in twTablePrimes(g);
    if not (trace_ok or weil_ok) or not twHasOps(X) then return true, _; end if;
    LD := twLevelData(X`D, X`N);
    // The ops are built on every prime either filter uses on X, not on p alone, so that the safety
    // net of twCurveData (each h must commute with T_l for every such l) is no weaker than theirs.
    L := X`D*X`N;
    ps := Sort(SetToSequence(Set(twTracePrimes(g, L)) join Set(twWeilPrimes(g, L)) join {p}));
    Tk, ops, IK := twCurveData(X, LD, ps);
    if trace_ok then
        ok, h, _, v, tr := twTraceCheckOn(X, Tk, ops, IK, [p]);
        if not ok then
            return false, Sprintf("TwistedTrace, h = %o with p^v = %o^%o", h, p, v);
        end if;
    end if;
    if weil_ok then
        TAB := twWeilTable(g, p);
        P1 := twWeilPolynomial(Tk[p], IK, IK, p, g);
        for o in ops do
            Ph := twWeilPolynomial(Tk[p], o[2], IK, p, g);
            assert ChangeRing(Ph, GF(2)) eq ChangeRing(P1, GF(2));
            if twWeilKey(Ph) notin TAB then
                return false, Sprintf("TwistedWeilPolynomial h = %o with p = %o", o[1], p);
            end if;
        end for;
    end if;
    return true, _;
end intrinsic;

intrinsic TwistedWeilPolynomials(X::ShimuraQuot, p::RngIntElt) -> SeqEnum
{The twisted Weil polynomials P_(X_h) at the good prime p, as a list of <name of h, P_h>, for h = 1
followed by the admissible involutions h != 1 of X defined over Q.  At h = 1 this is the Weil
polynomial of X.}
    require (X`D*X`N) mod p ne 0 : "p must be a good prime";
    LD := twLevelData(X`D, X`N);
    Tk, ops, IK := twCurveData(X, LD, [p]);
    return [<"1", twWeilPolynomial(Tk[p], IK, IK, p, X`g)>] cat
           [<o[1], twWeilPolynomial(Tk[p], o[2], IK, p, X`g)> : o in ops];
end intrinsic;

// The undecided curves of genus >= 3 in curves, grouped by level (D,N) in increasing D*N; only the
// curves with some prime to test (per primes_of(g, D*N)) and some admissible involution h != 1.
function twLevels(curves, primes_of)
    levels := AssociativeArray();
    for i->X in curves do
        if assigned X`IsSubhyp then continue; end if;
        if X`g lt 3 then continue; end if;
        if #primes_of(X`g, X`D*X`N) eq 0 then continue; end if;
        if not twHasOps(X) then continue; end if;
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
        LD := twLevelData(key[1], key[2]);
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
        LD := twLevelData(key[1], key[2]);
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
