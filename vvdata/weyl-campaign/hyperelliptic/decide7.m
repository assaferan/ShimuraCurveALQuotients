// The seven undecided curves whose bases have model files: build the fibre product of the stored
// double covers (coordinate check, genus, trace formula at three good primes), then Noether's count.
AttachSpec("ShimuraQuotients.spec"); SetColumns(0);
P<x> := PolynomialRing(Rationals());
curves := GetHyperellipticCandidates();
todo := [<390,1,[1,10,26,65]>, <390,1,[1,6,26,39]>, <462,1,[1,6,11,66]>, <510,1,[1,2,15,30]>, <510,1,[1,2,17,34]>, <22,15,[1,2,55,110]>, <10,51,[1,15,102,170]>];
for t in todo do
    D := t[1]; N := t[2]; W := {Integers()| w : w in t[3]};
    X := rep{Y : Y in curves | Y`D eq D and Y`N eq N and Y`W eq W};
    models := eval (Read(Sprintf("data/models/models_%o_%o.m", D, N)) cat "\nreturn models;");
    covers := AssociativeArray();
    DN := D*N; full := {Integers()| d : d in Divisors(DN) | GCD(d, DN div d) eq 1};
    for k in Keys(models) do
        if #models[k] gt 0 and Type(models[k][1][2]) ne MonStgElt and #Set(k)*2 eq #full then
            e := models[k][1]; covers[Set(k)] := e[2] + (#e ge 3 select e[3]^2/4 else 0);
        end if;
    end for;
    failing, classes := AtkinLehnerTripleCheck(covers, D, N);
    ups := AtkinLehnerDoubleCoversOver(W, D, N);
    ok, gens := FibreProductGenerators(W, {U : U in ups | IsDefined(covers, U)}, D, N);
    assert ok;
    bad, tri := DefinitelyDisagreeingFactors(gens, failing, classes);
    fs := [covers[U] : U in gens];
    K := FibreProductFunctionField(fs);
    g := Genus(K);
    // trace formula at three good primes
    nok := 0; p := 3; tried := 0;
    while nok lt 3 and tried lt 40 do
        p := NextPrime(p); tried +:= 1; if DN mod p eq 0 then continue; end if;
        Kp := RationalFunctionField(GF(p)); L := Kp;
        try for f in fs do R<Y> := PolynomialRing(L); L := FunctionField(Y^2 - L!Evaluate(PolynomialRing(GF(p))!f, Kp.1)); end for; catch err continue; end try;
        if Genus(L) ne X`g then continue; end if;
        cnt := [&+[ee * #Places(L, ee) : ee in Divisors(d)] : d in [1..2]];
        assert cnt eq [ComputePointsViaTrace(X, p, d) : d in [1..2]];
        nok +:= 1;
    end while;
    hyp, nrel := IsHyperellipticFunctionField(K);
    okc, c := ConicFactorOfFibreProduct(fs);
    printf "X_0(%o,%o)/%o: genus %o (expected %o), factors %o, triple check %o, trace formula at %o primes; %o (%o quadrics)%o\n",
        D, N, t[3], g, X`g, [Sort(SetToSequence(U)) : U in gens], bad select "DISAGREE" else TripleCheckVerdict(gens, failing, classes), nok,
        hyp select "HYPERELLIPTIC" else "NOT hyperelliptic", nrel,
        hyp select Sprintf("; conic factor %o, rational point %o", okc select c else "none", okc select ConicHasRationalPoint(c) else "-") else "";
end for;
quit;
