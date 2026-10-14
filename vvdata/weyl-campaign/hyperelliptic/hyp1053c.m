AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
P<x> := PolynomialRing(Rationals());
curves := GetHyperellipticCandidates();
models := eval (Read("data/models/models_10_53.m") cat "\nreturn models;");
for k in [[1,530], [1,53]] do
    e := models[[Integers()|k[1], k[2]]][1];
    fs := CRVDoubleCovers(e[3]);
    K := FibreProductFunctionField(fs);
    hyp, nrel := IsHyperellipticFunctionField(K);
    ok, c, S := ConicFactorOfFibreProduct(fs);
    printf "%o: genus %o hyperelliptic %o (%o quadrics); conic factor %o with rational point %o\n", k, Genus(K), hyp, nrel, c, ConicHasRationalPoint(c);
    H := HyperellipticModelOverConic(fs);
    f, h := HyperellipticPolynomials(H);
    X_ := rep{Y : Y in curves | Y`D eq 10 and Y`N eq 53 and Y`W eq {Integers()|k[1], k[2]}};
    nok := 0;
    for p in [3, 7, 11, 13, 17, 19, 23, 29] do
        if not IsNonsingular(ChangeRing(H, GF(p))) then continue; end if;
        assert [#Points(ChangeRing(H, GF(p^d))) : d in [1..2]] eq [ComputePointsViaTrace(X_, p, d) : d in [1..2]];
        nok +:= 1;
    end for;
    printf "   y^2 = %o\n   genus %o, trace formula at %o primes\n   MODEL %o\n", f, Genus(H), nok, Eltseq(f);
end for;
quit;
