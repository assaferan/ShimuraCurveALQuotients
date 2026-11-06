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
// Jacobian to be isotypic for a newform old at 2, which is rechecked.

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
        M, B, variant, md5 := Explode(spec);
        fn := Sprintf("%o/basis_%o_%o.txt", DualTestBasisDirectory(), M, B);
        if not FileExists(fn) then return false, {}, Sprintf("the weight-3/2 basis %o is missing", fn); end if;
        if variant ne "sharp" then
            ok, why := jacobian_old_at_two(star);
            if not ok then return false, {}, why; end if;
        end if;
        Append(~tests, DualTestSpace(D, N, ReadWeightThreeHalvesBasis(M, B : MD5 := md5) : Variant := variant));
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
    t := JacobianTorsionBound(star : MaxPrime := cert["torsion"][1], NumPrimes := cert["torsion"][2]);
    if #eff le t then return false, Sprintf("%o divisors against a torsion bound of %o", #eff, t); end if;
    return true, Sprintf("HyperellipticCertificate A: %o degree-two CM divisors related by %o, torsion bound %o (%o)",
                         #eff, info, t, cert["file"]);
end function;

intrinsic VerifyHyperellipticCertificate(X::ShimuraQuot, cert::Assoc) -> BoolElt, MonStgElt
{True and a label if the certificate cert, read from a file of data/certificates/, proves that X is
 hyperelliptic; otherwise false and the hypothesis that fails.  Every hypothesis is rechecked.}
    if <cert["D"], cert["N"]> ne <X`D, X`N> or {Integers() | m : m in cert["W"]} ne X`W then
        return false, "the certificate is for another curve";
    end if;
    case cert["shape"]:
        when "A": return verify_shape_A(X, cert);
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
{Marks hyperelliptic every undecided curve of genus at least three with a certificate in Dir that
 VerifyHyperellipticCertificate accepts.  A refused certificate is reported and the curve left
 undecided.  Prints one summary line.}
    nok := 0; nref := 0;
    for i->X in curves do
        if assigned X`IsSubhyp or X`g lt 3 then continue; end if;
        for cert in ReadHyperellipticCertificates(X`D, X`N : Dir := Dir) do
            if {Integers() | m : m in cert["W"]} ne X`W then continue; end if;
            ok, label := VerifyHyperellipticCertificate(X, cert);
            if ok then
                curves[i]`IsSubhyp := true;
                curves[i]`IsHyp := true;
                curves[i]`TestInWhichProved := label;
                nok +:= 1;
                break;
            end if;
            printf "HyperellipticCertificates: X_0^%o(%o)/%o: certificate refused: %o\n", X`D, X`N, Sort(SetToSequence(X`W)), label;
            nref +:= 1;
        end for;
    end for;
    printf "HyperellipticCertificates summary: %o curve(s) proved hyperelliptic, %o certificate(s) refused\n", nok, nref;
end intrinsic;
