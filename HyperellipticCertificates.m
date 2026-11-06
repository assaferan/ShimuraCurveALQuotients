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
// full Atkin-Lehner quotient.

// The star X_0^D(N)/W_full as the record the Borcherds routines expect.
function cert_star(D, N)
    full := {Integers() | m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1};
    X := CreateShimuraQuot(D, N, full);
    X`g := GenusShimuraCurveQuotient(D, N, full);
    X`CurveID := 0;
    _ := exists(star){c : c in GetQuotientsAndGenera([X]) | IsStarCurve(c)};
    return star, full;
end function;

function verify_shape_A(X, cert)
    star, full := cert_star(X`D, X`N);
    if X`W ne full then return false, "shape A needs the full Atkin-Lehner group"; end if;
    eff := cert["divisors"];
    // effective, of degree two, distinct
    for E in eff do
        if exists{c : c in E | c[2] le 0} then return false, Sprintf("divisor %o is not effective", E); end if;
        deg := &+[c[2] * Degree(FieldOfDefinitionOfCMPoint(star, c[1])) : c in E];
        if deg ne 2 then return false, Sprintf("divisor %o has degree %o", E, deg); end if;
    end for;
    sets := [{<c[1], c[2]> : c in E} : E in eff];
    if #Set(sets) ne #eff then return false, "two of the divisors are equal"; end if;
    supp := [{c[1] : c in E} : E in eff];
    // each form's divisor is a difference D_i - D_j of two of them, with disjoint supports
    fs := BorcherdsFormsFromFile(cert["forms"], star);
    edges := {};
    for k in Keys(fs) do
        dv := Set(DivisorOfBorcherdsForm(fs[k], star));
        hit := [<i, j> : i, j in [1..#eff] | i ne j and IsEmpty(supp[i] meet supp[j]) and
                dv eq sets[i] join {<c[1], -c[2]> : c in sets[j]}];
        if IsEmpty(hit) then return false, Sprintf("form %o has divisor %o, not a difference of the divisors", k, dv); end if;
        Include(~edges, {hit[1][1], hit[1][2]});
        // an integral multiple of the form is a Borcherds input: integral principal part, no constant term
        g := qExpansionAtoo(fs[k], 1);
        cf := [Coefficient(g, m) : m in [Valuation(g)..0]];
        n := LCM([Denominator(c) : c in cf]);
        if Coefficient(g, 0) ne 0 or not forall{c : c in cf | IsIntegral(n*c)} then
            return false, Sprintf("form %o: no integral multiple of its principal part", k);
        end if;
    end for;
    // the differences connect all the divisors, so every [D_i - D_1] is torsion
    reached := {1};
    repeat
        new := reached join {Rep(e diff reached) : e in edges | #(e meet reached) eq 1};
        grown := #new gt #reached; reached := new;
    until not grown;
    if #reached ne #eff then return false, "the differences given by the forms do not connect the divisors"; end if;
    t := JacobianTorsionBound(star : MaxPrime := cert["torsion"][1], NumPrimes := cert["torsion"][2]);
    if #eff le t then return false, Sprintf("%o divisors against a torsion bound of %o", #eff, t); end if;
    return true, Sprintf("HyperellipticCertificate A: %o degree-two CM divisors related by %o Borcherds forms, torsion bound %o (%o)",
                         #eff, #Keys(fs), t, cert["file"]);
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
