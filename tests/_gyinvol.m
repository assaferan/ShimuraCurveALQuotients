// tests/_gyinvol.m -- NOT a test (leading underscore: excluded from the CI matrix). Run it by hand
// to REGENERATE the `ws_data` matrices pasted into tests/X0_{55_1,22_3,15_2,14_5}.m. It is kept
// because it is the PROVENANCE of those matrices: without it they are unexplained integers.
//
// Transport Guo-Yang's PUBLISHED Atkin-Lehner involutions into OUR stored model's coordinates.
//
// NON-CIRCULAR: the involutions are external (Guo-Yang's table); the coordinate change is computed
// from the two EQUATIONS alone (IsIsomorphic), never from the pipeline's own `ws` output. The
// pipeline's involutions are what tests/BorcherdsProducts.m then checks against the result.
//
// The psi it picks is one element of a torsor under Aut; a different choice conjugates ALL the
// transported involutions simultaneously, and the harness searches that same torsor, so the
// ambiguity is absorbed rather than being a source of false verdicts.
AttachSpec("ShimuraQuotients.spec");

P<x> := PolynomialRing(Rationals());

// base, key, GY's f, [<label, 3x3 matrix in GY's coords, row-vector convention>]
data := [*
  <"55_1", [Integers()|1], -(x^4-x^3+x^2+x+1)*(3*x^4+x^3-5*x^2-x+3),
    [* <5,  Matrix(Rationals(),3,3,[0,0,1, 0,1,0, -1,0,0])>,      // (-1/x, y/x^4)
       <55, DiagonalMatrix(Rationals(),[1,-1,1])> *]>,            // (x, -y)
  <"22_3", [Integers()|1], -27*x^8 - 308*x^6 - 2146*x^4 - 308*x^2 - 27,
    [* <2,  Matrix(Rationals(),3,3,[0,0,1, 0,-1,0, -1,0,0])>,     // (-1/x, -y/x^4)
       <3,  DiagonalMatrix(Rationals(),[-1,1,1])>,                // (-x, y)
       <66, DiagonalMatrix(Rationals(),[1,-1,1])> *]>,            // (x, -y)
  <"15_2", [Integers()|1], -(x^2+3)*(3*x^2+4)*(x^4-x^2+4),
    [* <2,  Matrix(Rationals(),3,3,[0,0,1, 0,-4,0, 2,0,0])>,      // (2/x, -4y/x^4)
       <3,  DiagonalMatrix(Rationals(),[-1,1,1])>,                // (-x, y)
       <5,  DiagonalMatrix(Rationals(),[-1,-1,1])> *]>,           // (-x, -y)
  <"14_5", [Integers()|1], -23*x^8 - 180*x^7 - 358*x^6 - 168*x^5 - 677*x^4
                            + 168*x^3 - 358*x^2 + 180*x - 23,
    [* <2,  Matrix(Rationals(),3,3,[0,0,1, 0,1,0, -1,0,0])>,      // (-1/x, y/x^4)
       <14, DiagonalMatrix(Rationals(),[1,-1,1])>,                // (x, -y)
       <35, Matrix(Rationals(),3,3,[1,0,2, 0,-25,0, 2,0,-1])> *]> // ((x+2)/(2x-1), -25y/(2x-1)^4)
*];

for d in data do
    b, key, fgy, invs := Explode(d);
    printf "\n================ %o ================\n", b;
    models := eval (Read("data/models/models_" cat b cat ".m") cat "\nreturn models;");
    ok, e := IsDefined(models, key);
    if not ok then printf "  NO KEY %o\n", key; continue; end if;
    fours := e[1][2];
    Cours := HyperellipticCurve(fours);
    Cgy   := HyperellipticCurve(fgy);
    printf "  our genus %o, GY genus %o\n", Genus(Cours), Genus(Cgy);
    if Genus(Cours) ne Genus(Cgy) then printf "  GENUS MISMATCH\n"; continue; end if;
    isit, psi := IsIsomorphic(Cours, Cgy);
    if not isit then printf "  NOT ISOMORPHIC\n"; continue; end if;
    Rg := CoordinateRing(AmbientSpace(Cgy));
    cg := [Rg.i : i in [1..Rank(Rg)]];
    for iv in invs do
        Q, M := Explode(iv);
        wgy := map< Cgy -> Cgy | Eltseq(Vector(cg)*ChangeRing(M, Rg)) >;
        // sanity: it must actually be an involution of GY's curve
        printf "  w_%o: GY-side involution? ", Q;
        printf "%o\n", (wgy*wgy eq IdentityMap(Cgy));
        ourw := psi*wgy*Inverse(psi);
        de := DefiningEquations(ourw);
        Ro := CoordinateRing(AmbientSpace(Cours));
        n := Rank(Ro);
        lin := true;
        Mo := ZeroMatrix(Rationals(), n, n);
        for j in [1..n] do
            pj := Ro!de[j];
            for m in Terms(pj) do
                ex := Exponents(LeadingMonomial(m));
                if &+ex ne 1 then lin := false; break j; end if;
                i := Index(ex, 1);
                Mo[i][j] := LeadingCoefficient(m);
            end for;
        end for;
        if not lin then printf "    w_%o NOT LINEAR: %o\n", Q, de; continue; end if;
        // it must be an involution OF OUR CURVE, expressed by that matrix
        co := [Ro.i : i in [1..n]];
        chk := map< Cours -> Cours | Eltseq(Vector(co)*ChangeRing(Mo, Ro)) >;
        printf "    ws_data[{1}][%o] := Matrix(%o,%o,%o);   // matrix-form == transported: %o; involution: %o\n",
               Q, n, n, [Rationals()| Mo[i][j] : j in [1..n], i in [1..n]],
               chk eq ourw, (chk*chk eq IdentityMap(Cours));
    end for;
end for;
exit;
