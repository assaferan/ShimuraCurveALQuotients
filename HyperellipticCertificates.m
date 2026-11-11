// HyperellipticCertificates.m
//
// A pipeline stage that decides a curve from a certificate file instead of a search: each file
// data/certificates/<D>_<N>.m holds the inputs of a proof that some X_0^D(N)/W is hyperelliptic,
// and the stage rechecks every hypothesis of that proof.  The searches that found the inputs are
// not repeated and the stage never depends on them.
//
// Shape A (torsion counting, paper/torsion-counting-hyperelliptic.tex, Theorem 1.2).  Let X have
// genus at least two and rational torsion of order at most t in its Jacobian.  If D_1, ..., D_(t+1)
// are distinct effective divisors of degree two on X and each D_i - D_j along a connected set of
// pairs is, up to a positive multiple, the divisor of a Borcherds form with integral principal
// part, then every class [D_i - D_1] is rational torsion (Borcherds, Duke Math. J. 97 (1999),
// Theorem 4.1, with the correction in Duke Math. J. 105 (2000)), two of them coincide, and X has a
// function of degree two.  The forms are Borcherds products on the star, so shape A applies to the
// full Atkin-Lehner quotient.  The relations come either from banked Borcherds forms
// (cert["forms"]) or from the dual test of DualTest.m against saved weight-3/2 bases
// (cert["dualtests"], cert["relations"]); the 2 || N variants of the dual test also need the
// Jacobian to be isotypic for a newform old at 2, which is rechecked.  The torsion bound may be lowered
// by component groups (cert["phi"], paper/component-groups.tex): an l-part of the rational torsion that
// injects into Phi_p is trivial on a torsion CM divisor whose class in Phi_p is zero.

// The star X_0^D(N)/W_full as the record the Borcherds routines expect.
function cert_star(D, N)
    full := {Integers() | m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1};
    X := CreateShimuraQuot(D, N, full);
    X`g := GenusShimuraCurveQuotient(D, N, full);
    X`CurveID := 0;
    _ := exists(star){c : c in GetQuotientsAndGenera([X]) | IsStarCurve(c)};
    return star, full;
end function;

// The relations from banked Borcherds forms: each form's divisor is a difference D_i - D_j with
// disjoint supports, and an integral multiple of the form is a Borcherds input.
function relations_by_forms(star, cert, sets, supp, ndiv)
    fs := BorcherdsFormsFromFile(cert["forms"], star);
    edges := {};
    for k in Keys(fs) do
        dv := Set(DivisorOfBorcherdsForm(fs[k], star));
        hit := [<i, j> : i, j in [1..ndiv] | i ne j and IsEmpty(supp[i] meet supp[j]) and
                dv eq sets[i] join {<c[1], -c[2]> : c in sets[j]}];
        if IsEmpty(hit) then return false, {}, Sprintf("form %o has divisor %o, not a difference of the divisors", k, dv); end if;
        Include(~edges, {hit[1][1], hit[1][2]});
        g := qExpansionAtoo(fs[k], 1);
        cf := [Coefficient(g, m) : m in [Valuation(g)..0]];
        n := LCM([Denominator(c) : c in cf]);
        if Coefficient(g, 0) ne 0 or not forall{c : c in cf | IsIntegral(n*c)} then
            return false, {}, Sprintf("form %o: no integral multiple of its principal part", k);
        end if;
    end for;
    return true, edges, Sprintf("%o Borcherds forms", #Keys(fs));
end function;

// The Hecke polynomial prod (x - a_p(f)) of the Jacobian of X at a good prime p, from the trace formula.
function hecke_poly(X, p)
    g := X`g;
    R<x> := PolynomialRing(Rationals());
    s := [ p^i + 1 - ComputePointsViaTrace(X, p, i) : i in [1..g] ];
    e := [Rationals() | 1];
    for k in [1..g] do Append(~e, (&+[ (-1)^(i-1) * e[k-i+1] * s[i] : i in [1..k] ]) / k); end for;
    E := func< k | k le g select e[k+1] else p^(k-g) * e[2*g-k+1] >;
    P := &+[ (-1)^k * E(k) * x^(2*g-k) : k in [0..2*g] ];
    Mat := ZeroMatrix(Rationals(), g+1, 2*g+1);
    for k in [0..g] do for j in [0..k] do Mat[k+1, g + k - 2*j + 1] +:= Binomial(k, j) * p^j; end for; end for;
    ok, q := IsConsistent(Mat, Vector(Rationals(), [ Coefficient(P, i) : i in [0..2*g] ]));
    assert ok;
    return R![ q[k+1] : k in [0..g] ];
end function;

// The hypothesis of the 2 || N dual tests: the Jacobian of the star is isotypic for one newform orbit
// of dimension g (its Hecke polynomial at one of the first six good primes is irreducible), and the
// only orbit of level D*m, m | N, with the Hecke polynomials of the Jacobian at those six primes has
// odd level.
function jacobian_old_at_two(star)
    D := star`D; N := star`N; g := star`g;
    ps := [p : p in PrimesInInterval(3, 200) | (D*N) mod p ne 0][1..6];
    polys := [hecke_poly(star, p) : p in ps];
    if not exists{P : P in polys | IsIrreducible(P)} then
        return false, Sprintf("no Hecke polynomial at %o is irreducible, so the Jacobian is not shown isotypic", ps);
    end if;
    matches := [];
    for m in Divisors(N) do
        Snew := NewSubspace(CuspidalSubspace(ModularSymbols(D*m, 2, +1)));
        for A in NewformDecomposition(Snew) do
            if Dimension(A) eq g and [HeckePolynomial(A, p) : p in ps] eq polys then Append(~matches, D*m); end if;
        end for;
    end for;
    if #matches ne 1 then return false, Sprintf("%o newform orbits match the Jacobian", #matches); end if;
    if IsEven(matches[1]) then return false, Sprintf("the Jacobian is new at 2 (level %o)", matches[1]); end if;
    return true, Sprintf("Jacobian old at 2, level %o", matches[1]);
end function;

// The relations from the dual test: relation <i, j, k> says that dual test k certifies D_i - D_j.
function relations_by_dual_test(star, cert, eff)
    D := star`D; N := star`N;
    tests := [* *];
    for spec in cert["dualtests"] do
        M, B, variant, md5 := Explode(<spec[1], spec[2], spec[3], spec[4]>);
        fn := Sprintf("%o/basis_%o_%o.txt", DualTestBasisDirectory(), M, B);
        if not FileExists(fn) then return false, {}, Sprintf("the weight-3/2 basis %o is missing", fn); end if;
        if variant ne "sharp" then
            ok, why := jacobian_old_at_two(star);
            if not ok then return false, {}, why; end if;
        end if;
        T := DualTestSpace(D, N, ReadWeightThreeHalvesBasis(M, B : MD5 := md5) : Variant := variant);
        if #spec ge 5 then T := RestrictDualTestSpace(T, spec[5][1], spec[5][2]); end if;
        Append(~tests, T);
    end for;
    edges := {};
    for r in cert["relations"] do
        i, j, k := Explode(r);
        target := eff[i] cat [<c[1], -c[2]> : c in eff[j]];
        if not CertifiesTorsion(tests[k], target) then
            return false, {}, Sprintf("the dual test at level %o does not certify D_%o - D_%o", cert["dualtests"][k][1], i, j);
        end if;
        Include(~edges, {i, j});
    end for;
    return true, edges, Sprintf("the dual test at level %o", Join([Sprint(t[1]) : t in cert["dualtests"]], ", "));
end function;

function verify_shape_A(X, cert)
    star, full := cert_star(X`D, X`N);
    if X`W ne full then return false, "shape A needs the full Atkin-Lehner group"; end if;
    eff := cert["divisors"];
    // effective, of degree two, distinct
    for E in eff do
        if exists{c : c in E | c[2] le 0} then return false, Sprintf("divisor %o is not effective", E); end if;
        deg := &+[c[2] * DegreeOfFieldOfDefinitionOfCMPoint(star, c[1]) : c in E];
        if deg ne 2 then return false, Sprintf("divisor %o has degree %o", E, deg); end if;
    end for;
    sets := [{<c[1], c[2]> : c in E} : E in eff];
    if #Set(sets) ne #eff then return false, "two of the divisors are equal"; end if;
    supp := [{c[1] : c in E} : E in eff];
    if IsDefined(cert, "forms") then
        ok, edges, info := relations_by_forms(star, cert, sets, supp, #eff);
    else
        ok, edges, info := relations_by_dual_test(star, cert, eff);
    end if;
    if not ok then return false, info; end if;
    // the differences connect all the divisors, so every [D_i - D_1] is torsion
    reached := {1};
    repeat
        new := reached join {Rep(e diff reached) : e in edges | #(e meet reached) eq 1};
        grown := #new gt #reached; reached := new;
    until not grown;
    if #reached ne #eff then return false, "the certified differences do not connect the divisors"; end if;
    // Check := false: the bound needs the counts over F_p, ..., F_(p^g) only; the further extension
    // the check adds needs class numbers of discriminants near 4 p^(g+1), out of reach at p ~ 100
    t := JacobianTorsionBound(star : MaxPrime := cert["torsion"][1], NumPrimes := cert["torsion"][2], Check := false);
    // cert["phi"]: primes p | D whose component groups remove l-parts.  Every D_i - D_1 is torsion; an l
    // whose part injects into Phi_p, where D_i - D_1 has class zero, does not divide its order.  For the l
    // removed from every D_i - D_1, all the classes lie in the subgroup with trivial l-parts.
    tt := t; how := "";
    if IsDefined(cert, "phi") and #eff ge 2 then
        common := {l : l in PrimeDivisors(t)};
        for i in [2..#eff] do
            Bi := TorsionOrderBoundOfCMDivisor(X`D, X`N, full, eff[i] cat [<c[1], -c[2]> : c in eff[1]] :
                                               Bound := t, Primes := cert["phi"]);
            common meet:= {l : l in PrimeDivisors(t) | Bi mod l ne 0};
        end for;
        for l in common do while tt mod l eq 0 do tt div:= l; end while; end for;
        if tt ne t then how := Sprintf(", lowered to %o by the component groups at %o", tt, cert["phi"]); end if;
    end if;
    if #eff le tt then return false, Sprintf("%o divisors against a torsion bound of %o", #eff, tt); end if;
    return true, Sprintf("HyperellipticCertificate A: %o degree-two CM divisors related by %o, torsion bound %o%o (%o)",
                         #eff, info, t, how, cert["file"]);
end function;

// Shape B (genus four over a genus-two star, paper/hyperelliptic-candidates.tex, Proposition
// prop:normalform).  Y = X_0^D(N)/H of genus 4 is a double cover of the star C (genus 2) branched at
// the images of the fixed points of the involutions outside H; Y is hyperelliptic exactly when that
// branch divisor D_1 is a fibre of the hyperelliptic map of C, i.e. canonical.  The certificate relates
// D_1 to other effective degree-two divisors D_2, ... of C by torsion differences (Borcherds forms or
// the dual test, as in shape A); when the Jacobian of C is torsion-free, D_1 is linearly equivalent to
// a distinct effective divisor of degree two, so it is canonical.

// The common part of shapes B and B2: Y of genus four over the genus-two star, the first divisor its branch
// divisor, and the certified relations between the divisors.
function shape_B_relations(X, cert)
    star, full := cert_star(X`D, X`N);
    if #full ne 2 * #X`W or star`g ne 2 or X`g ne 4 then
        return false, "shape B needs a genus-four double cover of a genus-two star", _, _, _, _;
    end if;
    eff := cert["divisors"];
    // the branch divisor: the CM points of the discriminants carrying fixed points of the involutions
    // outside the group of Y
    branch := {Integers() | };
    for m in full diff X`W do
        e := NumFixedPointsByCMOrder(X`D, X`N, m);
        branch join:= {d : d in Keys(e) | e[d] ne 0};
    end for;
    if branch ne {c[1] : c in eff[1]} or exists{c : c in eff[1] | c[2] ne 1} then
        return false, Sprintf("the branch divisor is supported on %o, not the first divisor %o", branch, eff[1]), _, _, _, _;
    end if;
    for E in eff do
        if exists{c : c in E | c[2] le 0} then return false, Sprintf("divisor %o is not effective", E), _, _, _, _; end if;
        deg := &+[c[2] * DegreeOfFieldOfDefinitionOfCMPoint(star, c[1]) : c in E];
        if deg ne 2 then return false, Sprintf("divisor %o has degree %o", E, deg), _, _, _, _; end if;
    end for;
    sets := [{<c[1], c[2]> : c in E} : E in eff];
    if #Set(sets) ne #eff or #eff lt 2 then return false, "the divisors must be at least two and distinct", _, _, _, _; end if;
    supp := [{c[1] : c in E} : E in eff];
    if IsDefined(cert, "forms") then
        ok, edges, info := relations_by_forms(star, cert, sets, supp, #eff);
    else
        ok, edges, info := relations_by_dual_test(star, cert, eff);
    end if;
    if not ok then return false, info, _, _, _, _; end if;
    if not exists{e : e in edges | 1 in e} then return false, "no certified difference involves the branch divisor", _, _, _, _; end if;
    return true, "", star, eff, edges, info;
end function;

// A canonical branch divisor makes Y the double cover given by a theta characteristic L over Q, and Y is
// hyperelliptic only when L is odd.  Galois preserves the quadratic form of L on J[2], whose orthogonal
// group has order 72 when L is even and 120 when it is odd, so a Frobenius of order 5 on J[2] makes L odd.
// Frob_p acts on J[2] by the L-polynomial mod 2, x^4 + a x^3 + b x^2 + a x + 1, of order 5 iff a, b are odd.
// Returns the first such p, or 0.
function frobenius_of_order_five(star : MaxPrime := 200)
    for p in PrimesInInterval(3, MaxPrime) do
        if (star`D * star`N) mod p eq 0 then continue; end if;
        P := WeilPolynomial(star, p);
        if IsOdd(Integers()!Coefficient(P, 1)) and IsOdd(Integers()!Coefficient(P, 2)) then return p; end if;
    end for;
    return 0;
end function;

function verify_shape_B(X, cert)
    ok, why, star, eff, edges, info := shape_B_relations(X, cert);
    if not ok then return false, why; end if;
    t := JacobianTorsionBound(star : MaxPrime := cert["torsion"][1], NumPrimes := cert["torsion"][2], Check := false);
    if t ne 1 then return false, Sprintf("the Jacobian of the star has torsion bound %o, not 1", t); end if;
    p5 := frobenius_of_order_five(star);
    if p5 eq 0 then
        return false, "no Frobenius of order 5 on the 2-torsion of the star below 200, so the theta characteristic may be even";
    end if;
    return true, Sprintf("HyperellipticCertificate B: the branch divisor %o is canonical on the genus-two star, related to %o other divisor(s) by %o, torsion-free (%o); the theta characteristic is odd (Frobenius of order 5 on J[2] at p = %o)",
                         eff[1], #{Rep(e diff {1}) : e in edges | 1 in e}, info, cert["file"], p5);
end function;

// The Hecke polynomial x^2 - s x + n of a genus-two curve C over Q at a good odd prime p, from its
// L-polynomial 1 - s T + (2p + n) T^2 - p s T^3 + p^2 T^4.
function hecke_poly_genus2(C, p)
    L := LPolynomial(ChangeRing(C, GF(p)));
    R<x> := PolynomialRing(Rationals());
    s := -Coefficient(L, 1); n := Coefficient(L, 2) - 2*p;
    return R![n, -s, 1];
end function;

// Whether the endomorphism with tangent matrix T (in the basis x^(i-1) dx/y of the simplified model f)
// maps the period lattice of y^2 = f into itself, to 10^-40 at 80 digits.
function tangent_preserves_lattice(f, T)
    CC := ComplexField(80);
    Pi := BigPeriodMatrix(AnalyticJacobian(ChangeRing(f, CC)));
    cj := func< A | Matrix(CC, Nrows(A), Ncols(A), [ComplexConjugate(z) : z in Eltseq(A)]) >;
    TP := ChangeRing(T, CC) * Pi;
    R := VerticalJoin(Pi, cj(Pi))^-1 * VerticalJoin(TP, cj(TP));
    return Max([Abs(z - Round(Real(z))) : z in Eltseq(R)]) lt 10^-40;
end function;

// Shape B2 (paper/hyperelliptic-6-341.tex): as shape B, but the star's Jacobian J may have a rational point
// t of order two.  Relations joining t + 1 divisors to the branch divisor B put two of them in one class,
// which is then canonical, so B ~ K or B ~ K + t.  The second is excluded on models.  The Hecke algebra
// acts on J through the maximal order O_K of a real quadratic field, and every abelian surface isogenous to
// J with O_K in its endomorphisms is J or J/<t>: odd isogenies are ideal quotients (the reducibility and
// congruence gcds are powers of 2), and in the Bruhat-Tits tree at the prime over 2 every Galois-fixed
// vertex has at most one fixed neighbour, as the rational torsion of the class has order at most 2.  The
// certificate names two curves (cert["stars"]) whose Jacobians are two such surfaces, distinct by their
// invariants; their real multiplication is CHIMP's certified endomorphism (cert["rm"], the tangent
// matrix, which is rechecked to preserve the period lattice) and their conductor at 2 CHIMP's proved
// exponent (cert["conductor2"]), from which modularity identifies their newform with J's.  Each has
// O_K with totally positive units squares, so one principal polarisation, and J with its polarisation is
// the Jacobian of one of them; on both, K + t is not twice a rational class, so no double cover over Q is
// branched along it.
function verify_shape_B2(X, cert)
    ok, why, star, eff, edges, info := shape_B_relations(X, cert);
    if not ok then return false, why; end if;
    D := X`D; N := X`N; DN := D*N;
    t := JacobianTorsionBound(star : MaxPrime := cert["torsion"][1], NumPrimes := cert["torsion"][2], Check := false);
    if t gt 2 then return false, Sprintf("the Jacobian of the star has torsion bound %o, more than 2", t); end if;
    comp := {1};
    repeat n0 := #comp; comp join:= &join[e : e in edges | not IsEmpty(e meet comp)]; until #comp eq n0;
    if #comp lt t + 1 then
        return false, Sprintf("the branch divisor is related to %o divisor(s), fewer than the torsion bound %o", #comp - 1, t);
    end if;
    if t eq 1 then return verify_shape_B(X, cert); end if;
    // the Hecke order: T_p0 generates the maximal order of a real quadratic field K of class number one with
    // a unit of norm -1, in which 2 ramifies
    p0 := cert["heckeprime"];
    hp := hecke_poly(star, p0);
    dK := Integers()!Discriminant(hp);
    if not (dK gt 0 and IsFundamentalDiscriminant(dK) and IsEven(dK)) then
        return false, Sprintf("T_%o has discriminant %o: not the maximal order of a real quadratic field ramified at 2", p0, dK);
    end if;
    K := QuadraticField(dK);
    if ClassNumber(K) ne 1 or Norm(FundamentalUnit(MaximalOrder(K))) ne -1 then
        return false, Sprintf("Q(sqrt %o) has class number %o or no unit of norm -1", dK, ClassNumber(K));
    end if;
    // Jacquet-Langlands: J is isogenous to the newforms of level DN new at D with signs -1 at D, +1 at N
    if not IsSquarefree(DN) then return false, "shape B2 needs squarefree DN"; end if;
    S := CuspidalSubspace(ModularSymbols(DN, 2, 1));
    V := NewSubspace(S, D);
    for p in PrimeDivisors(DN) do V := AtkinLehnerSubspace(V, p, D mod p eq 0 select -1 else 1); end for;
    dec := [G : G in NewformDecomposition(NewSubspace(S)) | Dimension(G) eq 2];
    if Dimension(V) ne 2 or not exists(f){G : G in dec | G eq V} then
        return false, "the Atkin-Lehner eigenspace of the star is not one newform orbit of dimension two";
    end if;
    cq := func< h | [Rationals() | c : c in Eltseq(h)] >;
    good := [p : p in PrimesInInterval(5, 60) | DN mod p ne 0];
    hv := [cq(CharacteristicPolynomial(HeckeOperator(f, p))) : p in good];
    if [cq(hecke_poly(star, p)) : p in good] ne hv then
        return false, "the trace formula of the star and the newform disagree";
    end if;
    // the curves: conductor, real multiplication, the newform, the Mordell-Weil group, and K + t
    PQ<x> := PolynomialRing(Rationals());
    invs := [];
    for i -> c in cert["stars"] do
        C := HyperellipticCurve(PQ!c[1], PQ!c[2]);
        f2 := HyperellipticPolynomials(SimplifiedModel(C));
        cond := Conductor(C);
        if exists{p : p in PrimeDivisors(cond) | p ne 2 and Valuation(cond, p) ne 2 * Valuation(DN, p)}
           or exists{p : p in PrimeDivisors(DN) | p ne 2 and Valuation(cond, p) ne 2}
           or cert["conductor2"][i] ne 2 * Valuation(DN, 2) then
            return false, Sprintf("curve %o does not have conductor (DN)^2", i);
        end if;
        T := Matrix(Rationals(), 2, 2, cert["rm"][i]);
        if T^2 ne (dK div 4) * IdentityMatrix(Rationals(), 2) or not tangent_preserves_lattice(f2, T) then
            return false, Sprintf("the endomorphism cited for curve %o is not sqrt %o on its period lattice", i, dK div 4);
        end if;
        hc := [cq(hecke_poly_genus2(C, p)) : p in good];
        if hc ne hv or exists{G : G in dec | G ne f and [cq(CharacteristicPolynomial(HeckeOperator(G, p))) : p in good] eq hc} then
            return false, Sprintf("the Hecke polynomials of curve %o do not single out the newform of the star", i);
        end if;
        J := Jacobian(SimplifiedModel(C));
        G, mw, f1, f3 := MordellWeilGroupGenus2(J);
        if not (f1 and f3) then return false, Sprintf("the Mordell-Weil group of curve %o is not proved", i); end if;
        tors := [G.j : j in [1..Ngens(G)] | Order(G.j) ne 0];
        if &*[Integers() | Order(g) : g in tors] gt 2 then return false, Sprintf("curve %o has rational torsion of order more than 2", i); end if;
        if #tors eq 1 then
            pts := Points(SimplifiedModel(C) : Bound := 1000);
            if #pts eq 0 then return false, Sprintf("no rational point found on curve %o", i); end if;
            P0 := pts[1]; iP0 := SimplifiedModel(C)![P0[1], -P0[2], P0[3]];
            // L^2 = O(B) with B ~ K + t and L = O(P0) + s needs 2s = [iota P0 - P0] + t
            v := (J![iP0, P0] + mw(tors[1])) @@ mw;
            if v in sub<G | [2*G.j : j in [1..Ngens(G)]]> then
                return false, Sprintf("on curve %o the class K + t is twice a rational class", i);
            end if;
        end if;
        Append(~invs, G2Invariants(C));
    end for;
    if #invs ne 2 or invs[1] eq invs[2] then return false, "the certificate must name two non-isomorphic curves"; end if;
    // odd isogenies: reducibility at p = 1 mod DN and the congruence between the two embeddings
    C1 := HyperellipticCurve(PQ!cert["stars"][1][1], PQ!cert["stars"][1][2]);
    hs := [<p, hecke_poly_genus2(C1, p)> : p in cert["reducibility"] | IsPrime(p) and p mod DN eq 1];
    red := GCD([Integers() | Evaluate(h[2], 1 + h[1]) : h in hs]);
    cong := GCD([Integers() | Discriminant(hecke_poly_genus2(C1, p)) : p in PrimesInInterval(5, 400) | DN mod p ne 0]);
    if #hs lt 2 or red eq 0 or exists{q : q in PrimeDivisors(red * cong) | q ne 2} then
        return false, Sprintf("the reducibility gcd %o or the congruence gcd %o has an odd prime factor", red, cong);
    end if;
    return true, Sprintf("HyperellipticCertificate B2: the branch divisor %o is canonical on the genus-two star, related to %o other divisor(s) by %o; the star is one of two curves, on both of which the other class allowed by the torsion admits no double cover (%o)",
                         eff[1], #comp - 1, info, cert["file"]);
end function;

// E exactly: the isogeny class from the traces up to the Sturm bound of level DN, the curve from its
// reduction types I_n at the primes cert["tamagawa"] (default: those of D), n being the order of the
// component group of the genus-one quotient Q = X_0^D(N)/WQ (the star, or an intermediate quotient).
// Also returns the component-group data at those primes and cert["phi"].
function identify_elliptic_quotient(D, N, Q, WQ, cert)
    E := MinimalModel(EllipticCurve(cert["E"]));
    if (D*N) mod Conductor(E) ne 0 then return false, _, _, Sprintf("%o has conductor %o, not dividing %o", cert["E"], Conductor(E), D*N); end if;
    sturm := Ceiling((D*N) * &*[1 + 1/q : q in PrimeDivisors(D*N)] / 6);
    for q in [q : q in PrimesInInterval(2, sturm) | (D*N) mod q ne 0] do
        if q + 1 - ComputePointsViaTrace(Q, q, 1) ne TraceOfFrobenius(E, q) then
            return false, _, _, Sprintf("the quotient and %o differ at %o", cert["E"], q);
        end if;
    end for;
    tam := IsDefined(cert, "tamagawa") select cert["tamagawa"] else PrimeDivisors(D);
    phis := AssociativeArray();
    for q in Set(tam) join Set(cert["phi"]) do phis[q] := ComponentGroupData(D, N, WQ, q); end for;
    // the dual graph gives the geometric component group, of order n for reduction of type I_n, n the
    // valuation of the minimal discriminant (the Tamagawa number agrees only when the reduction is split,
    // as it is on the star, whose group contains w_p)
    orders := [&*([1] cat ComponentGroup(phis[q])) : q in tam];
    nI := func< F, q | Valuation(Discriminant(MinimalModel(F)), q) >;
    if [nI(E, q) : q in tam] ne orders then
        return false, _, _, Sprintf("the reduction types I_n of %o at %o are not the orders %o of the component groups", cert["E"], tam, orders);
    end if;
    if exists{E2 : E2 in IsogenousCurves(E) | not IsIsomorphic(E2, E) and [nI(E2, q) : q in tam] eq orders} then
        return false, _, _, "the component groups do not single out one curve of the isogeny class";
    end if;
    return true, E, phis, "";
end function;

function identify_star_curve(D, N, star, full, cert)
    return identify_elliptic_quotient(D, N, star, full, cert);
end function;

// Whether the rational point P of the minimal model E reduces to the identity component at p: to O when p
// divides a denominator, else to a smooth point of the reduction.
function in_identity_component(E, P, p)
    if IsZero(P) or exists{c : c in Eltseq(P) | Valuation(Denominator(c), p) gt 0} then return true; end if;
    Ep := ChangeRing(E, GF(p));
    return not (Curve(Ep)!(Ep(GF(p))![GF(p)!c : c in Eltseq(P)]) in SingularPoints(Curve(Ep)));
end function;

// Shape C (genus three over an elliptic star, paper/hyperelliptic-15-38.tex, Sections 1-4).  Y = X_0^D(N)/H
// of genus 3 is a double cover of the genus-one star E branched at four points; Y is hyperelliptic when they
// split into two pairs A, B of equal sum in the group E and Y is the untwisted cover.  The certificate names
// E (a Cremona label), the halves A and B, the dual test that makes A - B torsion, the primes p | D whose
// component groups Phi_p kill the torsion, and, when E(Q)[2] is not zero, a prime at which the Prym of
// Y -> E has an irreducible L-polynomial (the twist by the rational 2-torsion class has a Prym isogenous to
// a product of two elliptic curves over Q, paper/hyperelliptic-6-385.tex).
function verify_shape_C(X, cert)
    D := X`D; N := X`N;
    star, full := cert_star(D, N);
    if #full ne 2 * #X`W or star`g ne 1 or X`g ne 3 then
        return false, "shape C needs a genus-three double cover of a genus-one star";
    end if;
    A := cert["divisors"][1]; B := cert["divisors"][2];
    // the branch divisor is A + B
    branch := {Integers() | };
    for m in full diff X`W do
        e := NumFixedPointsByCMOrder(D, N, m);
        branch join:= {d : d in Keys(e) | e[d] ne 0};
    end for;
    if branch ne {c[1] : c in A cat B} or exists{c : c in A cat B | c[2] ne 1} or #(A cat B) ne #branch then
        return false, Sprintf("the branch divisor is supported on %o, not on the halves %o, %o", branch, A, B);
    end if;
    for E0 in [A, B] do
        deg := &+[c[2] * DegreeOfFieldOfDefinitionOfCMPoint(star, c[1]) : c in E0];
        if deg ne 2 then return false, Sprintf("the half %o has degree %o", E0, deg); end if;
    end for;
    ok, E, phis, why := identify_star_curve(D, N, star, full, cert);
    if not ok then return false, why; end if;
    // A - B is torsion.  On a genus-one star every Heegner component lies in the a_q(E)-eigenspace of
    // T(q^2), so a test spec <M, B, variant, md5, q> is restricted to it, with a_q taken from E.
    eff := [A, B];
    cc := cert;
    cc["dualtests"] := [ #spec ge 5 select <spec[1], spec[2], spec[3], spec[4], <spec[5], TraceOfFrobenius(E, spec[5])>>
                                    else spec : spec in cert["dualtests"] ];
    ok, edges, info := relations_by_dual_test(star, cc, eff);
    if not ok then return false, info; end if;
    if {1, 2} notin edges then return false, "the dual test does not relate the two halves"; end if;
    // A - B = 0: its class is zero in each Phi_p of cert["phi"], and every nonzero rational torsion point
    // leaves the identity component at one of those p
    Delta := A cat [<c[1], -c[2]> : c in B];
    for q in cert["phi"] do
        if exists{x : x in ClassInComponentGroup(phis[q], Delta) | x ne 0} then
            return false, Sprintf("the pairing has nonzero class in Phi_%o", q);
        end if;
    end for;
    G, mG := TorsionSubgroup(E);
    if exists{T : T in G | T ne G!0 and forall{q : q in cert["phi"] | in_identity_component(E, mG(T), q)}} then
        return false, Sprintf("a rational torsion point of %o lies on the identity component at every p in %o", cert["E"], cert["phi"]);
    end if;
    // the twist
    two := #TwoTorsionSubgroup(E);
    how := "no rational 2-torsion";
    if two gt 1 then
        if IsDefined(cert, "prym") then
            q := cert["prym"];
            P := ZetaNumeratorViaTrace(X, q : Check := false) div ZetaNumeratorViaTrace(star, q : Check := false);
            if not IsIrreducible(P) then return false, Sprintf("the Prym L-polynomial at %o is reducible", q); end if;
            how := Sprintf("the Prym simple (L-polynomial irreducible at %o)", q);
        elif IsDefined(cert, "dualgraph") then
            if two ne 2 then return false, "the dual-graph twist test needs E(Q)[2] = Z/2"; end if;
            q, P0 := Explode(cert["dualgraph"]);
            if not IsDefined(phis, q) then phis[q] := ComponentGroupData(D, N, full, q); end if;
            ok, why := CoverIsUntwistedAtPrime(phis[q], X`W, A, B, P0);
            if not ok then return false, why; end if;
            how := why;
        else
            return false, "E(Q)[2] is not zero and no twist argument is given";
        end if;
    end if;
    return true, Sprintf("HyperellipticCertificate C: the halves %o and %o have equal sum on %o, torsion by %o, zero by Phi_p for p in %o; untwisted: %o (%o)",
                         A, B, cert["E"], info, cert["phi"], how, cert["file"]);
end function;

// Shape C4 (genus three over an elliptic star, branched at one place of degree four; paper/hyperelliptic-6-323.tex).
// The four points of discriminant d are a torsor under Cl(d)/Cl(d)[2] = Z/4, the only Galois-stable pairing is
// {P, s^2 P} | {s P, s^3 P}, and Delta = P + s^2 P - s P - s^3 P is the CM cycle twisted by the genus character
// chi = (d1, d2), d1 d2 = d.  It is torsion by the Gross-Zagier formula of Yuan-Zhang-Zhang when
// L(E/Q(sqrt d), chi, s) = L(E^(d1), s) L(E^(d2), s) vanishes to order at least three at s = 1: one twist has
// root number -1, the other root number +1 and L(1) = 0 exactly (the twisted winding element of E).  Delta lies in
// the minus part of E(K)_tors for K one of Q(sqrt d1), Q(sqrt d2), which is {0, T}; at p | D the four points lie
// on two components, two each, so the class of Delta in Phi_p over K_p is 0 or plus or minus twice their
// difference, and when that is neither 0 nor the image of T, Delta = 0.  The cover is untwisted by a simple Prym.
function verify_shape_C4(X, cert)
    D := X`D; N := X`N;
    star, full := cert_star(D, N);
    if #full ne 2 * #X`W or star`g ne 1 or X`g ne 3 then
        return false, "shape C4 needs a genus-three double cover of a genus-one star";
    end if;
    d := cert["d"];
    branch := {Integers() | };
    for m in full diff X`W do
        e := NumFixedPointsByCMOrder(D, N, m);
        branch join:= {dd : dd in Keys(e) | e[dd] ne 0};
    end for;
    if branch ne {d} then return false, Sprintf("the branch divisor is supported on %o, not on %o", branch, d); end if;
    if DegreeOfFieldOfDefinitionOfCMPoint(star, d) ne 4 then return false, "the branch place does not have degree four"; end if;
    ok, E, phis, why := identify_star_curve(D, N, star, full, cert);
    if not ok then return false, why; end if;
    // the Galois structure and the genus character
    Cl, mCl := ClassGroup(QuadraticForms(d));
    amb := sub< Cl | [ g : g in Cl | 2*g eq Cl!0 ] >;
    Q := quo< Cl | amb >;
    if not (IsCyclic(Q) and #Q eq 4) then return false, "Cl(d)/Cl(d)[2] is not cyclic of order four"; end if;
    d1, d2 := Explode(cert["chi"]);
    if d1 * d2 ne d or not IsFundamentalDiscriminant(d1) or not IsFundamentalDiscriminant(d2) then
        return false, Sprintf("(%o, %o) is not a genus character of discriminant %o", d1, d2, d);
    end if;
    chi := func< q | d1 mod q eq 0 select KroneckerSymbol(d2, q) else KroneckerSymbol(d1, q) >;
    if exists{q : q in PrimeDivisors(d) | chi(q) ne 1} then return false, "the genus character is not trivial on the ambiguous classes"; end if;
    // the L-functions.  L(E^(dd), 1) is the twisted period sum_a chi_dd(a) <f, {0, a/|dd|}> up to a nonzero
    // constant (Birch's formula for the twist of f by a primitive character; at the primes of dd that
    // divide the level, E is multiplicative and both sides have no p-term), read off the modular symbols
    // of E at its own level, so the zero test is exact
    function twisted_L_vanishes(dd)
        eps := KroneckerCharacter(dd, Rationals());
        A := ModularSymbols(E, Evaluate(eps, -1) eq 1 select 1 else -1);
        return IsZero(RationalMapping(A)(TwistedWindingElement(AmbientSpace(A), 1, eps)));
    end function;
    w := []; vanish := [];
    for dd in [d1, d2] do
        Et := MinimalModel(QuadraticTwist(E, dd));
        Append(~w, RootNumber(Et));
        Append(~vanish, RootNumber(Et) eq -1 select 1 else (twisted_L_vanishes(dd) select 2 else 0));
    end for;
    if &+vanish lt 3 then return false, Sprintf("L(E/k, chi, s) is not shown to vanish to order 3 (root numbers %o)", w); end if;
    // the minus part of the torsion over both candidate fields
    G, mG := TorsionSubgroup(E);
    two := [ mG(g) : g in G | g ne G!0 and 2*g eq G!0 ];
    if #two ne 1 then return false, "E(Q)[2] is not Z/2"; end if;
    T := two[1];
    for dd in [d1, d2] do
        K := QuadraticField(dd); EK := BaseChange(E, K);
        GK, mK := TorsionSubgroup(EK);
        minus := [ g : g in GK | g ne GK!0 and EK![ Conjugate(c) : c in Eltseq(mK(g)) ] eq -mK(g) ];
        if #minus ne 1 or mK(minus[1]) ne EK!T then return false, Sprintf("the minus part of E(Q(sqrt %o))_tors is not {0, T}", dd); end if;
    end for;
    // the component groups at p
    p := cert["phi"][1];
    if in_identity_component(E, T, p) then return false, Sprintf("T reduces to the identity component at %o", p); end if;
    Gp := phis[p];
    cnt, total := CMSpecialisation(Gp, d);
    vs := Sort(Setseq(Keys(cnt)));
    if total ne 4 or #vs ne 2 or exists{v : v in vs | cnt[v] ne 2 or v gt Gp`nv} then
        return false, "the four points do not lie two each on two components";
    end if;
    function critical(th)
        V := Gp`nv; Ed := [];
        for i->e in Gp`edges do
            prev := e[1];
            for k in [1..th[i]-1] do V +:= 1; Append(~Ed, <prev, V>); prev := V; end for;
            Append(~Ed, <prev, e[2]>);
        end for;
        Lap := ZeroMatrix(Integers(), V, V);
        for e in Ed do Lap[e[1],e[1]] +:= 1; Lap[e[2],e[2]] +:= 1; Lap[e[1],e[2]] -:= 1; Lap[e[2],e[1]] -:= 1; end for;
        Sm, Pm := SmithForm(Lap);
        idx := [i : i in [1..V] | Sm[i,i] notin {0, 1}];
        if #idx ne 1 then return false, 0, 0; end if;
        n := Sm[idx[1], idx[1]];
        v := ZeroMatrix(Integers(), V, 1); v[vs[1], 1] := 1; v[vs[2], 1] := -1;
        return true, n, (Pm * v)[idx[1], 1] mod n;
    end function;
    for dd in [d1, d2] do
        K := QuadraticField(dd); OK := Integers(K);
        fac := Factorization(p * OK);
        if #fac eq 2 then
            ok, n, delta := critical(Gp`thick);
        elif fac[1][2] eq 2 then
            ok, n, delta := critical([2 * t : t in Gp`thick]);
            if ok and LocalInformation(EK, fac[1][1])[4] ne n where EK := BaseChange(E, K) then
                return false, "the component group over the ramified extension is not the Tamagawa number there";
            end if;
        else
            return false, Sprintf("%o is inert in Q(sqrt %o)", p, dd);
        end if;
        if not ok then return false, "the component group is not cyclic"; end if;
        if (2 * delta) mod n in {0, n div 2} then
            return false, Sprintf("over Q_%o(sqrt %o), twice the difference of the components is 0 or the image of T", p, dd);
        end if;
    end for;
    // untwisted
    q := cert["prym"];
    P := ZetaNumeratorViaTrace(X, q : Check := false) div ZetaNumeratorViaTrace(star, q : Check := false);
    if not IsIrreducible(P) then return false, Sprintf("the Prym L-polynomial at %o is reducible", q); end if;
    return true, Sprintf("HyperellipticCertificate C4: the place of discriminant %o on %o; Delta torsion by Yuan-Zhang-Zhang (root numbers %o, L(E^(%o),1) = 0 by the twisted modular symbol), zero by Phi_%o; untwisted: the Prym simple at %o (%o)",
                         d, cert["E"], w, w[1] eq 1 select d1 else d2, p, q, cert["file"]);
end function;

// Shape N (NOT hyperelliptic; a genus-three double cover of an elliptic star).  If Y were hyperelliptic, its
// hyperelliptic involution, defined over Q and commuting with the deck involution, would descend to x -> c - x
// on E and pair the four branch points Galois-stably into two pairs of equal sum (paper/hyperelliptic-15-38.tex,
// Section 1).  The verifier lists the Galois-stable pairings from the degrees of the branch places and shows
// that each difference A - B has nonzero class in a component group Phi_p, so is not zero in E.
function verify_shape_N(X, cert)
    D := X`D; N := X`N;
    star, full := cert_star(D, N);
    if #full ne 2 * #X`W or star`g ne 1 or X`g ne 3 then
        return false, "shape N needs a genus-three double cover of a genus-one star";
    end if;
    branch := {Integers() | };
    for m in full diff X`W do
        e := NumFixedPointsByCMOrder(D, N, m);
        branch join:= {d : d in Keys(e) | e[d] ne 0};
    end for;
    places := Sort(Setseq(branch));
    degs := [DegreeOfFieldOfDefinitionOfCMPoint(star, d) : d in places];
    if &+degs ne 4 then return false, Sprintf("the branch places %o have degrees %o, not of total 4", places, degs); end if;
    rat := [places[i] : i in [1..#places] | degs[i] eq 1];
    quad := [places[i] : i in [1..#places] | degs[i] eq 2];
    // the Galois-stable pairings, as differences A - B
    pairings := [];
    if Sort(degs) eq [1, 1, 2] then
        pairings := [ [<rat[1], 1>, <rat[2], 1>, <quad[1], -1>] ];
    elif degs eq [1, 1, 1, 1] then
        r := rat;
        pairings := [ [<r[1], 1>, <r[j], 1>] cat [<r[k], -1> : k in [2..4] | k ne j] : j in [2..4] ];
    elif Sort(degs) eq [1, 3] then
        return true, Sprintf("HyperellipticCertificate N: not hyperelliptic, the branch places %o of degrees %o admit no Galois-stable pairing (%o)", places, degs, cert["file"]);
    else
        return false, Sprintf("branch places of degrees %o are not handled", degs);
    end if;
    used := [];
    for Delta in pairings do
        hit := 0;
        for p in cert["phi"] do
            G := ComponentGroupData(D, N, full, p);
            try
                c := ClassInComponentGroup(G, Delta);
            catch err
                continue;
            end try;
            if exists{x : x in c | x ne 0} then hit := p; break; end if;
        end for;
        if hit eq 0 then return false, Sprintf("the pairing %o has class zero in every Phi_p for p in %o", Delta, cert["phi"]); end if;
        Append(~used, hit);
    end for;
    return true, Sprintf("HyperellipticCertificate N: not hyperelliptic, every Galois-stable pairing of the branch places %o has nonzero class in Phi_p (p = %o) (%o)",
                         places, used, cert["file"]);
end function;

// Shape T (NOT hyperelliptic; a genus-three double cover of an elliptic intermediate quotient, decided by the
// twist; paper/classification-certificates.tex, Section 7).  E = X/<W, w_m> (m = cert["over"]) has genus one and
// Y -> E is branched at two places Z_1, Z_2 of degree two with different fields of definition, so Z_1 | Z_2 is the
// only Galois-stable pairing.  Each Z_i is all the CM points of its discriminant on E, two of them, and a further
// Atkin-Lehner involution t of E with genus-zero quotient moves them (a point of X has stabiliser of order at most
// two in the Atkin-Lehner group, and theirs is in W), so Z_i is a fibre over a rational point of the quotient and
// Z_1 - Z_2 is the divisor of a function x.  If Y were hyperelliptic it would be E(sqrt(c x)); the dual graph of Y
// at the odd prime p = cert["p"] dividing D shows it is not.
function verify_shape_T(X, cert)
    D := X`D; N := X`N; L := D*N; m := cert["over"]; p := cert["p"];
    full := {Integers() | q : q in Divisors(L) | GCD(q, L div q) eq 1};
    if not (IsPrime(p) and p ge 3 and D mod p eq 0) then return false, Sprintf("%o is not an odd prime dividing D", p); end if;
    if m notin full or m in X`W then return false, Sprintf("w_%o is not an Atkin-Lehner involution outside W", m); end if;
    WE := X`W join {AtkinLehnerMul(m, w, L) : w in X`W};
    E := CreateShimuraQuot(D, N, WE); E`g := GenusShimuraCurveQuotient(D, N, WE); E`CurveID := 0;
    if X`g ne 3 or E`g ne 1 then return false, "shape T needs a genus-three double cover of a genus-one quotient"; end if;
    branch := {Integers() | };
    for q in WE diff X`W do
        e := NumFixedPointsByCMOrder(D, N, q);
        branch join:= {d : d in Keys(e) | e[d] ne 0};
    end for;
    places := Sort(Setseq(branch));
    if #places ne 2 or exists{d : d in places | DegreeOfFieldOfDefinitionOfCMPoint(E, d) ne 2} then
        return false, Sprintf("the branch places %o are not two places of degree two", places);
    end if;
    F1 := FieldOfDefinitionOfCMPoint(E, places[1]); F2 := FieldOfDefinitionOfCMPoint(E, places[2]);
    if IsIsomorphic(F1, F2) then return false, "the two branch places have the same field of definition"; end if;
    G := ComponentGroupData(D, N, WE, p);
    for d in places do
        _, tot := CMSpecialisation(G, d);
        if tot ne 2 then return false, Sprintf("E has %o CM points of discriminant %o, not two", tot, d); end if;
    end for;
    ts := [t : t in Sort(Setseq(full diff WE)) | GenusShimuraCurveQuotient(D, N, WE join {AtkinLehnerMul(t, w, L) : w in WE}) eq 0];
    if #ts eq 0 then return false, "no Atkin-Lehner involution of E has a quotient of genus zero"; end if;
    t := ts[1];
    ok, why := CoverIsNotSquareRootAtPrime(G, X`W, [<places[1], 1>], [<places[2], 1>]);
    if not ok then return false, why; end if;
    return true, Sprintf("HyperellipticCertificate T: not hyperelliptic; over E = X/<W, w_%o> the only Galois-stable pairing of the branch places %o has equal sums (fibres of E/w_%o), and Y is not the square-root cover: %o (%o)",
                         m, places, t, why, cert["file"]);
end function;

// Shape V (NOT hyperelliptic; a Klein-four cover of a genus-zero star; paper/classification-certificates.tex,
// Section 7).  W has index four in the Atkin-Lehner group, the star C has genus zero and the three intermediate
// quotients genus one.  Over E = X/<W, w_m> (m = cert["over"]) the branch points are two places Z_1, Z_2 of
// degree two, of different discriminants and with different fields, so Z_1 | Z_2 is the only Galois-stable
// pairing.  Each Z_i is a fibre of E -> C (a point of X has stabiliser of order at most two in the Atkin-Lehner
// group), so a hyperelliptic Y would be E(sqrt(c x)) with x pulled back from C, and K(C)(sqrt(c x)) would be a
// quadratic subfield of K(Y)/K(C) of genus zero; the three quadratic subfields of that extension have genus one.
function verify_shape_V(X, cert)
    D := X`D; N := X`N; L := D*N; m := cert["over"];
    if not IsSquarefree(L) then return false, "shape V needs D*N squarefree"; end if;
    full := {Integers() | q : q in Divisors(L) | GCD(q, L div q) eq 1};
    if #full ne 4 * #X`W then return false, "W does not have index four in the Atkin-Lehner group"; end if;
    if m notin full or m in X`W then return false, Sprintf("w_%o is not an Atkin-Lehner involution outside W", m); end if;
    if X`g ne 3 or GenusShimuraCurveQuotient(D, N, full) ne 0 then
        return false, "shape V needs a curve of genus three over a star of genus zero";
    end if;
    mids := {X`W join {AtkinLehnerMul(t, w, L) : w in X`W} : t in full diff X`W};
    gs := [GenusShimuraCurveQuotient(D, N, Wt) : Wt in mids];
    if #mids ne 3 or exists{g : g in gs | g ne 1} then
        return false, Sprintf("the intermediate quotients have genera %o, not 1, 1, 1", gs);
    end if;
    WE := X`W join {AtkinLehnerMul(m, w, L) : w in X`W};
    E := CreateShimuraQuot(D, N, WE); E`g := 1; E`CurveID := 0;
    // a point of X_0^D(N) fixed by some w_q, q in WE \ W, has trivial stabiliser in W, so |W| of them lie over
    // one branch point of Y -> E
    cnt := AssociativeArray();
    for q in WE diff X`W do
        e := NumFixedPointsByCMOrder(D, N, q);
        for d in [d : d in Keys(e) | e[d] ne 0] do cnt[d] := (IsDefined(cnt, d) select cnt[d] else 0) + e[d]; end for;
    end for;
    places := Sort(Setseq(Keys(cnt)));
    if #places ne 2 or exists{d : d in places | cnt[d] ne 2 * #X`W or DegreeOfFieldOfDefinitionOfCMPoint(E, d) ne 2} then
        return false, Sprintf("the branch points over E, of discriminants %o, are not two places of degree two", places);
    end if;
    if IsIsomorphic(FieldOfDefinitionOfCMPoint(E, places[1]), FieldOfDefinitionOfCMPoint(E, places[2])) then
        return false, "the two branch places have the same field of definition";
    end if;
    return true, Sprintf("HyperellipticCertificate V: not hyperelliptic; Y is a Klein-four cover of a genus-zero star with three elliptic intermediate quotients, and over E = X/<W, w_%o> its branch places %o have different fields, so they are fibres over the star and a function of degree two on Y would give the star a genus-zero double cover inside Y (%o)",
                         m, places, cert["file"]);
end function;

// Shape O (a genus-three double cover of an elliptic star whose form is old at a prime p || N, p not dividing
// D; paper/classification-certificates.tex, Section 6).  The four branch points Q_1, ..., Q_4 are one Galois
// orbit of CM points of fundamental discriminant d, permuted by a D4 whose unique stable pairing
// {Q_1, Q_3} | {Q_2, Q_4} is defined over the quadratic subfield k, and R - R^sigma, R = Q_1 + Q_3, is the
// chi-component of the orbit for the genus character chi of k.  The star maps to E by phi o pi_1 + phi o pi_2,
// and pi_2 of a CM point is pi_1 of it moved by the Frobenius of a prime above p (p splits in Q(sqrt d)), so
// that component is multiplied by 1 + chi(p) = 0: R - R^sigma is torsion in E(k)^- = E^(d_k)(Q), which is
// torsion-free for every curve of the class, so the pairing has equal sums, and E(Q)[2] = 0 leaves no twist.
function verify_shape_O(X, cert)
    D := X`D; N := X`N; p := cert["p"];
    star, full := cert_star(D, N);
    if #full ne 2 * #X`W or star`g ne 1 or X`g ne 3 then
        return false, "shape O needs a genus-three double cover of a genus-one star";
    end if;
    if not (IsPrime(p) and N mod p eq 0 and (N div p) mod p ne 0 and D mod p ne 0) then
        return false, Sprintf("%o is not a prime exactly dividing N and prime to D", p);
    end if;
    // the star's Jacobian is the curve of conductor DN/p: its form is old at p
    E := MinimalModel(EllipticCurve(cert["E"]));
    if Conductor(E) ne (D*N) div p then return false, Sprintf("%o does not have conductor %o", cert["E"], (D*N) div p); end if;
    sturm := Ceiling((D*N) * &*[1 + 1/q : q in PrimeDivisors(D*N)] / 6);
    for q in [q : q in PrimesInInterval(2, sturm) | (D*N) mod q ne 0] do
        if q + 1 - ComputePointsViaTrace(star, q, 1) ne TraceOfFrobenius(E, q) then
            return false, Sprintf("the star and %o differ at %o", cert["E"], q);
        end if;
    end for;
    // the branch points: one place of degree 4, of a fundamental discriminant d
    branch := {Integers() | };
    for m in full diff X`W do
        e := NumFixedPointsByCMOrder(D, N, m);
        branch join:= {dd : dd in Keys(e) | e[dd] ne 0};
    end for;
    if #branch ne 1 then return false, Sprintf("the branch points have the discriminants %o, not one", branch); end if;
    d := rep{dd : dd in branch};
    if not IsFundamentalDiscriminant(d) or DegreeOfFieldOfDefinitionOfCMPoint(star, d) ne 4 then
        return false, Sprintf("the branch place of discriminant %o is not of degree 4 with d fundamental", d);
    end if;
    L := AbsoluteField(FieldOfDefinitionOfCMPoint(star, d));
    if #GaloisGroup(L) ne 8 then return false, "the Galois group of the branch place is not D4"; end if;
    ks := Subfields(L, 2);
    if #ks ne 1 then return false, "the field of the branch place does not have one quadratic subfield"; end if;
    dk := Discriminant(MaximalOrder(ks[1][1]));
    if dk eq d or d mod dk ne 0 or not IsFundamentalDiscriminant(d div dk) then
        return false, Sprintf("Q(sqrt %o) is not the field of a genus character of discriminant %o", dk, d);
    end if;
    if KroneckerSymbol(d, p) ne 1 or KroneckerSymbol(dk, p) ne -1 then
        return false, Sprintf("chi(p) is not -1: p = %o, d = %o, d_k = %o", p, d, dk);
    end if;
    for F in IsogenousCurves(E) do
        if #TwoTorsionSubgroup(F) ne 1 then return false, Sprintf("%o has rational 2-torsion", CremonaReference(F)); end if;
        if #TorsionSubgroup(QuadraticTwist(F, dk)) ne 1 then
            return false, Sprintf("the twist of %o by %o has torsion", CremonaReference(F), dk);
        end if;
    end for;
    // a rational point on the star, so a Galois-invariant divisor class is the class of a rational divisor
    d0 := cert["ratpoint"];
    if not (CarriesCMPoint(D, N, d0) and DegreeOfFieldOfDefinitionOfCMPoint(star, d0) eq 1) then
        return false, Sprintf("the star has no rational CM point of discriminant %o", d0);
    end if;
    return true, Sprintf("HyperellipticCertificate O: the branch place of discriminant %o over a D4 field with quadratic subfield Q(sqrt %o); the form of the star is old at %o, where chi = %o, so the stable pairing has equal sums on %o; untwisted: no rational 2-torsion (%o)",
                         d, dk, p, KroneckerSymbol(dk, p), cert["E"], cert["file"]);
end function;

// Shape D (the exact Heegner-height route, paper/heegner-height-models.tex with the certified mode of
// ModelByHeegnerHeights): Y of genus 3 over the elliptic star E, identified exactly as in shape C; the
// branch points' heights from Cai-Shu-Tian with certified L-values, every recognition proved, give a finite
// list of candidate branch divisors containing the true one; every double cover over each is built exactly
// and all but the ones consistent with the trace formula are excluded.  Y is one of the survivors, so when
// every survivor is hyperelliptic (Noether's count 5) Y is, and when none is, Y is not.
function verify_shape_D(X, cert)
    D := X`D; N := X`N;
    star, full := cert_star(D, N);
    // E = X_0^D(N)/Wp: the star by default, or the intermediate quotient cert["Wp"]
    Wp := IsDefined(cert, "Wp") select {Integers() | m : m in cert["Wp"]} else full;
    Q := star;
    if Wp ne full then
        Q := CreateShimuraQuot(D, N, Wp); Q`g := GenusShimuraCurveQuotient(D, N, Wp); Q`CurveID := 0;
    end if;
    if not (X`W subset Wp and #Wp eq 2 * #X`W and Q`g eq 1 and X`g eq 3) then
        return false, "shape D needs a genus-three double cover of a genus-one quotient";
    end if;
    ok, E, phis, why := identify_elliptic_quotient(D, N, Q, Wp, cert);
    if not ok then return false, why; end if;
    res := ModelByHeegnerHeights(X, Wp, E : Certified := true, Effort := IsDefined(cert, "effort") select cert["effort"] else 1);
    if #res eq 0 then return false, "no configuration of branch points is consistent with the trace formula"; end if;
    noe := {r[5] : r in res};
    if noe eq {5} then
        return true, Sprintf("HyperellipticCertificate D: on %o, the certified Heegner heights leave %o configuration(s) consistent with the trace formula, all hyperelliptic (%o)", cert["E"], #res, cert["file"]);
    end if;
    return false, Sprintf("the surviving configurations have Noether counts %o", noe);
end function;

intrinsic VerifyHyperellipticCertificate(X::ShimuraQuot, cert::Assoc) -> BoolElt, MonStgElt
{True and a label if the certificate cert, read from a file of data/certificates/, proves its claim about X:
 that X is hyperelliptic, or for shapes N, T and V that it is not; otherwise false and the hypothesis that fails.
 Every hypothesis is rechecked.}
    if <cert["D"], cert["N"]> ne <X`D, X`N> or {Integers() | m : m in cert["W"]} ne X`W then
        return false, "the certificate is for another curve";
    end if;
    case cert["shape"]:
        when "A": return verify_shape_A(X, cert);
        when "B": return verify_shape_B(X, cert);
        when "B2": return verify_shape_B2(X, cert);
        when "C": return verify_shape_C(X, cert);
        when "C4": return verify_shape_C4(X, cert);
        when "N": return verify_shape_N(X, cert);
        when "T": return verify_shape_T(X, cert);
        when "V": return verify_shape_V(X, cert);
        when "D": return verify_shape_D(X, cert);
        when "O": return verify_shape_O(X, cert);
    end case;
    return false, Sprintf("shape %o is not implemented", cert["shape"]);
end intrinsic;

intrinsic ReadHyperellipticCertificates(D::RngIntElt, N::RngIntElt : Dir := "data/certificates") -> SeqEnum
{The certificates stored for X_0^D(N) in Dir/<D>_<N>.m, each with its file name under "file"; empty
 when there is no file.}
    fn := Sprintf("%o/%o_%o.m", Dir, D, N);
    if not FileExists(fn) then return []; end if;
    certs := eval (Read(fn) cat "\nreturn certificates;");
    for i in [1..#certs] do certs[i]["file"] := fn; end for;
    return certs;
end intrinsic;

intrinsic FilterByHyperellipticCertificates(~curves::SeqEnum : Dir := "data/certificates")
{Decides every undecided curve of genus at least three with a certificate in Dir that
 VerifyHyperellipticCertificate accepts: hyperelliptic, or not hyperelliptic for a certificate of shape N, T or V.  A refused certificate is reported and the curve left
 undecided.  Prints one summary line.}
    nok := 0; nnot := 0; nref := 0;
    for i->X in curves do
        if assigned X`IsSubhyp or X`g lt 3 then continue; end if;
        for cert in ReadHyperellipticCertificates(X`D, X`N : Dir := Dir) do
            if {Integers() | m : m in cert["W"]} ne X`W then continue; end if;
            ok, label := VerifyHyperellipticCertificate(X, cert);
            if ok then
                hyp := cert["shape"] notin {"N", "T", "V"};
                curves[i]`IsSubhyp := hyp;
                curves[i]`IsHyp := hyp;
                curves[i]`TestInWhichProved := label;
                if hyp then nok +:= 1; else nnot +:= 1; end if;
                break;
            end if;
            printf "HyperellipticCertificates: X_0^%o(%o)/%o: certificate refused: %o\n", X`D, X`N, Sort(SetToSequence(X`W)), label;
            nref +:= 1;
        end for;
    end for;
    printf "HyperellipticCertificates summary: %o curve(s) proved hyperelliptic, %o proved not hyperelliptic, %o certificate(s) refused\n", nok, nnot, nref;
end intrinsic;
