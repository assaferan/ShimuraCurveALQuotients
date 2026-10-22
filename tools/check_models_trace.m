// tools/check_models_trace.m
//
// Every stored model of positive genus in one model file against the Eichler-Selberg trace formula
// at every prime p with PMIN <= p <= PMAX not dividing DN: the numbers of points over F_p and
// F_{p^2} of the stored curve must equal ComputePointsViaTrace.  Hyperelliptic entries y^2 + h y = f
// are counted on the curve; fibre-product entries ("CRV") through the function field of the tower
// of double covers y_i^2 = f_i(s).  A prime where the model has bad reduction (a singular reduction
// or p in a denominator) is skipped and counted as such.
//
// Point counts at a few primes cannot tell a model from its quadratic twist by d when the character
// of d is -1 exactly at the primes where the trace vanishes; enough primes separate every twist by
// a d supported on the primes of DN and -1.
//
//   magma -b D:=33 N:=10 MODELS:=data/models/models_33_10.m PMIN:=13 PMAX:=113 tools/check_models_trace.m < /dev/null
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
P<x> := PolynomialRing(Rationals());
D := StringToInteger(D); N := StringToInteger(N);
pmin := assigned PMIN select StringToInteger(PMIN) else 13;
pmax := assigned PMAX select StringToInteger(PMAX) else 113;
models := eval (Read(MODELS) cat "\nreturn models;");
curves := GetHyperellipticCandidates();
primes := [p : p in PrimesInInterval(pmin, pmax) | (D*N) mod p ne 0];

// the double covers y_i^2 = f_i(s) of a stored fibre-product entry: each string is y_i^2 + F_i(s, z)
function covers_from_strings(strs)
    n := #strs;
    R := PolynomialRing(Rationals(), 2 + n);
    s := R.1; z := R.2;
    ys := [R.(2+i) : i in [1..n]];
    fs := [];
    for str in strs do
        F := eval (&cat[Sprintf("%o := R.%o; ", nm, j) : j->nm in ["s", "z"] cat [Sprintf("y%o", i) : i in [1..n]]] cat "return " cat str cat ";");
        i := [j : j in [1..n] | Degree(F, ys[j]) gt 0];
        assert #i eq 1 and F - ys[i[1]]^2 eq Evaluate(F, [s, z] cat [R!0 : j in [1..n]]);
        G := Evaluate(F - ys[i[1]]^2, [P.1, P!1] cat [P!0 : j in [1..n]]);
        Append(~fs, -G);                                    // y^2 + G = 0  <=>  y^2 = -G
    end for;
    return fs;
end function;

nok := 0; nbad := 0; nskip := 0; ncurves := 0;
for k in Sort([k : k in Keys(models)]) do
    if #models[k] eq 0 then continue; end if;
    e := models[k][1];
    g := e[1];
    if g eq 0 then continue; end if;
    X := rep{Y : Y in curves | Y`D eq D and Y`N eq N and Y`W eq Set(k)};
    assert X`g eq g;
    ncurves +:= 1;
    crv := Type(e[2]) eq MonStgElt;
    if crv then fs := covers_from_strings(e[3]); end if;
    for p in primes do
        cnt := [];
        try
            if crv then
                Kp := RationalFunctionField(GF(p)); L := Kp;
                for f in fs do
                    R<Y> := PolynomialRing(L);
                    L := FunctionField(Y^2 - L!Evaluate(PolynomialRing(GF(p))!f, Kp.1));
                end for;
                if Genus(L) ne g then nskip +:= 1; continue; end if;
                cnt := [&+[ee * #Places(L, ee) : ee in Divisors(d)] : d in [1..2]];
            else
                f := e[2]; h := #e ge 3 select e[3] else P!0;
                C := HyperellipticCurve(PolynomialRing(GF(p))!f, PolynomialRing(GF(p))!h);
                if Genus(C) ne g then nskip +:= 1; continue; end if;
                cnt := [#C, #BaseChange(C, GF(p^2))];
            end if;
        catch err
            nskip +:= 1; continue;
        end try;
        exp := [ComputePointsViaTrace(X, p, d) : d in [1..2]];
        if cnt eq exp then nok +:= 1; else nbad +:= 1; printf "MISMATCH %o_%o %o p=%o model %o trace %o\n", D, N, k, p, cnt, exp; end if;
    end for;
end for;
printf "%o: %o curves of positive genus, %o primes (%o..%o): %o prime-counts agree (F_p and F_p^2 each), %o disagree, %o skipped (bad reduction)\n",
    MODELS, ncurves, #primes, pmin, pmax, nok, nbad, nskip;
quit;
