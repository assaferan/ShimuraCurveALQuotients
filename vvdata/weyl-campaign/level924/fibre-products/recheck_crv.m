// Re-read the fibre-product entries of X_0^22(21) and X_0^6(77) from their STORED STRINGS alone and
// compare genus and point counts over F_p, F_{p^2} with the trace formula at the 4th and 5th good
// primes -- primes the writer (which used the first three) never looked at.
//   magma -b recheck_crv.m < /dev/null     (from a tree holding the appended model files)
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
P<x> := PolynomialRing(Rationals());
curves := GetHyperellipticCandidates();
nok := 0; nbad := 0; nskip := 0;
for b in [<22, 21>, <6, 77>] do
    D := b[1]; N := b[2]; DN := D*N;
    models := eval (Read(Sprintf("data/models/models_%o_%o.m", D, N)) cat "\nreturn models;");
    for k in Keys(models) do
        if #models[k] eq 0 or Type(models[k][1][2]) ne MonStgElt then continue; end if;
        e := models[k][1];
        g := e[1]; strs := e[3];
        nf := #strs;
        // the fibre-product entries have one y_i per factor and no x; the older pairs use x, y, s, z
        if exists{t : t in strs | "x" in t and not "y" in t} or exists{t : t in strs | " y " in t or "y^" in t or "y*" in t} then nskip +:= 1; continue; end if;
        assert nf le 6;
        R := PolynomialRing(Rationals(), 8);
        AssignNames(~R, ["s", "z"] cat ["y" cat IntegerToString(i) : i in [1..6]]);
        // eval sees the enclosing scope, so the coordinate names must exist as variables here
        s := R.1; z := R.2; y1 := R.3; y2 := R.4; y3 := R.5; y4 := R.6; y5 := R.7; y6 := R.8;
        fs := [];
        ok := true;
        for i in [1..nf] do
            E := eval strs[i];
            F := R.(2+i)^2 - E;             // the stored equation is y_i^2 - F_i(s, z)
            if exists{j : j in [3..8] | Degree(F, R.j) gt 0} then ok := false; break; end if;
            Append(~fs, P!Evaluate(F, [x, 1, 0, 0, 0, 0, 0, 0]));
        end for;
        if not ok then printf "SKIP W=%o: a stored equation is not of the form y_i^2 - F(s,z)\n", k; nskip +:= 1; continue; end if;
        X := rep{Y : Y in curves | Y`D eq D and Y`N eq N and Y`W eq Set(k)};
        assert X`g eq g;
        good := [p : p in PrimesUpTo(200) | DN mod p ne 0];
        checked := 0; idx := 3;      // start after the first three good primes
        while checked lt 2 and idx lt #good do
            idx +:= 1; p := good[idx];
            Kp := RationalFunctionField(GF(p)); L := Kp;
            try
                for f in fs do
                    fp := PolynomialRing(GF(p))!f;
                    RY<Y> := PolynomialRing(L);
                    L := FunctionField(Y^2 - L!Evaluate(fp, Kp.1));
                end for;
            catch err
                continue;
            end try;
            if Genus(L) ne g then continue; end if;
            cnt := [&+[d * #Places(L, d) : d in Divisors(m)] : m in [1..2]];
            exp := [ComputePointsViaTrace(X, p, m) : m in [1..2]];
            if cnt eq exp then nok +:= 1; else nbad +:= 1; printf "MISMATCH X_0^%o(%o) W=%o genus %o p=%o: counts %o, trace formula %o\n", D, N, k, g, p, cnt, exp; end if;
            checked +:= 1;
        end while;
        if checked lt 2 then printf "X_0^%o(%o) W=%o: only %o usable primes among the 4th..%oth good ones\n", D, N, k, checked, idx; end if;
    end for;
end for;
printf "fibre-product entries re-read from stored strings: %o (prime, cover) point counts agree with the trace formula, %o disagree; %o non-fibre-product CRV entries skipped\n", nok, nbad, nskip;
