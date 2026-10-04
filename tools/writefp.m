// tools/writefp.m -- the writer that produced the fibre-product entries of 2026-10-03.
//
// For every cover X_0(D,N)/W with no committed entry whose index-2 Atkin-Lehner double covers
// y^2 + h y = f(t) are committed in the same model file (over the star), it builds the compositum
// of the double covers (FibreProductCovers.m), keeps it only if the genus equals the Shimura-curve
// genus formula AND the point counts over F_p and F_{p^2} at the first three good primes agree with
// the Eichler-Selberg trace formula (a prime where some factor is a square is skipped, so an entry
// may rest on two primes; every generating set of double covers is tried, degenerate ones skipped),
// and writes the entry as a line in the existing non-hyperelliptic shape
// <genus, "CRV", [ Strings() | eq, ... ]> to OUTDIR/models_D_N.append, one file per base.  Ambient
// weights are not stored: s, z have weight 1 and a fibre coordinate has half the degree of its own
// equation, as the committed x, y, s, z pairs do.  Every rejection is printed with its reason.
//
//     NORMALIZ_BIN=... magma -b OUTDIR:=/some/dir tools/writefp.m < /dev/null
//
// ⚠ The factors are taken from the committed file by key; a file that mixes equations from runs
// with different Hauptmodul normalisations (the rebase stage can leave one so) gives a compositum
// of the right genus and the wrong curve, which only the trace-formula check rejects.
// ⚠ At level N > 1 the committed double covers depend on the m = 0 term of Schofer's formula
// (PR #66, a proposed proof), so these entries are checked by point counts, not proved.
AttachSpec("ShimuraQuotients.spec");
_ := ClassNumberLU(-4);
SetColumns(0);
P<x> := PolynomialRing(Rationals());
curves := eval Read("data/curves_after_UpdateCurves8.dat");
OUT := OUTDIR;
has := AssociativeArray(); ent := AssociativeArray();
files := Split(Pipe("ls data/models/models_*.m", ""), "\n");
for f in files do
    if #f eq 0 then continue; end if;
    parts := Split(Split(f, "/")[#Split(f, "/")], "_");
    D := StringToInteger(parts[2]); N := StringToInteger(Split(parts[3], ".")[1]);
    models := eval (Read(f) cat "\nreturn models;");
    s := {}; anyk := {};
    for k in Keys(models) do
        if #models[k] gt 0 then Include(~anyk, Set(k)); end if;
        if #models[k] gt 0 and Type(models[k][1][2]) ne MonStgElt then
            Include(~s, Set(k));
            // y^2 + h y = f is (y + h/2)^2 = f + h^2/4
            ent[<D,N,Set(k)>] := models[k][1][2] + (#models[k][1] ge 3 select models[k][1][3]^2/4 else 0);
        end if;
    end for;
    has[<D,N>] := <s, anyk>;
end for;
written := AssociativeArray();
nbuilt := 0; nrej := 0;
for X in curves do
    key := <X`D, X`N>;
    if not IsDefined(has, key) then continue; end if;
    DN := X`D * X`N;
    full := {Integers()| d : d in Divisors(DN) | GCD(d, DN div d) eq 1};
    if X`W eq full or X`W in has[key][2] then continue; end if;   // star, or already has SOME entry
    avail := {U : U in AtkinLehnerDoubleCoversOver(X`W, X`D, X`N) | U in has[key][1]};
    // every generating set of double covers is tried, as the stage does (FibreProductCovers.m); a
    // set some of whose subsets multiply to a constant times a square is degenerate and skipped
    kgen := Ilog2(#full div #X`W);
    if #avail lt kgen then continue; end if;
    genlist := [* *];
    for c in Subsets(avail, kgen) do
        I := full; for W2 in c do I := I meet W2; end for;
        if I eq X`W then Append(~genlist, SetToSequence(c)); end if;
    end for;
    if IsEmpty(genlist) then continue; end if;
    found := false; reasons := [];
    for gens in genlist do
    fs := [ent[<X`D, X`N, U>] : U in gens];
    if IsDegenerateFactorSet(fs) then
        Append(~reasons, "degenerate factor set"); continue;
    end if;
    K := FibreProductFunctionField(fs);
    if Genus(K) ne X`g then
        Append(~reasons, Sprintf("genus %o, expected %o", Genus(K), X`g)); continue;
    end if;
    // trace formula at three GOOD primes, F_p and F_{p^2}: a prime where a factor is a square mod p
    // or the genus drops is bad reduction of this model and is replaced by the next prime
    agree := true; checked := 0; p := 3; tried := 0;
    while checked lt 3 and tried lt 60 do
        p := NextPrime(p); tried +:= 1;
        if DN mod p eq 0 then continue; end if;
        Kp := RationalFunctionField(GF(p)); L := Kp;
        try
            for f in fs do
                R<Y> := PolynomialRing(L);
                L := FunctionField(Y^2 - L!Evaluate(PolynomialRing(GF(p))!f, Kp.1));
            end for;
        catch err
            continue;
        end try;
        if Genus(L) ne X`g then continue; end if;
        cnt := [&+[e * #Places(L, e) : e in Divisors(d)] : d in [1..2]];
        exp := [ComputePointsViaTrace(X, p, d) : d in [1..2]];
        if cnt ne exp then
            Append(~reasons, Sprintf("right genus, trace formula disagrees at p = %o", p));
            agree := false; break;
        end if;
        checked +:= 1;
    end while;
    if not agree then continue; end if;
    if checked lt 3 then Append(~reasons, Sprintf("only %o good primes among the first 60", checked)); continue; end if;
    found := true; break;
    end for;
    if not found then
        printf "REJECT %o_%o W=%o: %o\n", X`D, X`N, Sort(SetToSequence(X`W)), reasons;
        nrej +:= 1; continue;
    end if;
    C, eqns := FibreProductCurve(fs);
    Wseq := Sort(SetToSequence(X`W));
    line := Sprintf("models[[Integers()|%o]] := [* <%o, \"CRV\", [ Strings() | %o ]> *];",
        &cat[Sprintf("%o%o", Wseq[i], i lt #Wseq select "," else "") : i in [1..#Wseq]],
        X`g, &cat[Sprintf("\"%o\"%o", eqns[i], i lt #eqns select ", " else "") : i in [1..#eqns]]);
    fn := Sprintf("%o/models_%o_%o.append", OUT, X`D, X`N);
    if not IsDefined(written, key) then
        written[key] := 0;
        Write(fn, Sprintf("// Built 2026-10-03 as fibre products over the star line of committed double covers\n// (FibreProductCovers.m); each verified against the trace formula at up to three primes, F_p and F_{p^2}.\n// Coordinates s, z of weight 1; y_i of weight half the degree of its own equation.") : Overwrite := true);
    end if;
    Write(fn, line);
    written[key] +:= 1; nbuilt +:= 1;
end for;
printf "wrote %o entries across %o bases; rejected %o\n", nbuilt, #Keys(written), nrej;
