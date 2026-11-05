// tests/LMFDB138.m
//
// X_0(6,23) AND ITS ATKIN-LEHNER QUOTIENTS AGAINST THE NEWFORMS OF LEVEL 138, FROM LMFDB.
//
// Every other check of these models compares them with the repo's own evaluation of the
// Eichler-Selberg trace formula.  This one takes the same curves to a source outside the repo.
// By Jacquet-Langlands the Jacobian of X_0(6,23) is isogenous to the 6-new part of J_0(138), which
// for 138 = 2*3*23 is the whole new part: the four weight-2 newforms of level 138 and trivial
// character, of dimensions 1, 1, 1 and 2, five in total, the genus of X_0(6,23).  The Jacobian of
// a quotient X_0(6,23)/W is the part of that space on which every w in W acts trivially, so
//
//     #(X_0(6,23)/W)(F_p) = p + 1 - sum over the selected f of Tr a_p(f).
//
// THE SIGN RULE.  For w = prod of primes p, the eigenvalue of w on f is the product of sigma_p(f),
// where sigma_23(f) is the Atkin-Lehner eigenvalue of f at 23 and sigma_p(f) is MINUS the
// Atkin-Lehner eigenvalue at p = 2 and p = 3, the primes of the discriminant D = 6: the
// Jacquet-Langlands correspondence reverses the sign of the Atkin-Lehner involution at a prime
// where the quaternion algebra ramifies.  The rule is not assumed here, it is what the test
// confirms -- with it every curve below matches at every prime, and the dimensions match the
// genera; with the unflipped sign almost nothing matches.
//
// PROVENANCE OF THE NUMBERS.  LMFDB newforms 138.2.a.a, .b, .c (dimension 1) and 138.2.a.d
// (dimension 2), fetched 2026-10-13 from https://www.lmfdb.org/api/mf_newforms/?level=138&weight=2
// (fields atkin_lehner_eigenvals and traces, the latter being Tr a_n for n = 1, ..., 32).  The
// transcription is checked below two ways before it is used: the traces must be multiplicative in
// the Hecke sense (which determines all 32 of them from the values at the 11 primes), and at the
// primes dividing the level a_p must be -dim * w_p.
//
// F_{p^2} IS CHECKED WHERE THE DATA REACH.  #C(F_{p^2}) = p^2 + 1 - sum over the selected f of
// (Tr a_{p^2}(f) - p * dim f).  For a form of dimension 1, a_{p^2} = a_p^2 - p; for the
// dimension-2 form the trace of a_{p^2} is not determined by the trace of a_p, and the fetched
// traces reach n = 32, so it is available only at p = 5.  Curves whose selected forms are all of
// dimension 1 are therefore checked over F_{p^2} at every prime, the others only at p = 5.
//
// WHICH CURVES.  Every entry of data/models/models_6_23.m of genus at least 1: the ten
// hyperelliptic ones y^2 + h y = f, the pair stored at W = {1,6}, and the two fibre products
// (W = {1} and W = {1,23}), read back from their stored equations alone; and the top curve
// X_0(6,23) once more, rebuilt from three of the stored double covers as the fibre product over
// the star line, the way tests/FibreProductCovers.m builds it.  Point counts are numbers of places
// of degree 1 (and of degree 1 and 2) of the function field, so they are counts on the smooth
// curve.  A prime at which the genus drops is bad reduction of the model and is skipped.
//
// WHAT MAKES IT FAIL.  With the Atkin-Lehner sign at 2 and 3
// left unflipped, the very first curve already has the wrong dimension (1 against 2), so the
// Jacquet-Langlands sign is what the agreement below is evidence for.  Replacing every stored
// hyperelliptic f by -f gives 27 disagreements in the point counts.

printf "Checking X_0(6,23) and its quotients against the LMFDB newforms of level 138...\n";

lmf_label := ["138.2.a.a", "138.2.a.b", "138.2.a.c", "138.2.a.d"];
lmf_dim   := [1, 1, 1, 2];
// Atkin-Lehner eigenvalues at 2, 3, 23
lmf_al    := [[1, 1, 1], [1, -1, 1], [-1, 1, 1], [-1, -1, -1]];
// Tr a_n for n = 1, ..., 32
lmf_tr    := [
 [1, -1, -1, 1, -2, 1, -2, -1, 1, 2, -6, -1, -2, 2, 2, 1, 0, -1, 0, -2, 2, 6, -1, 1, -1, 2, -1, -2, 6, -2, 8, -1],
 [1, -1, 1, 1, 0, -1, 2, -1, 1, 0, 0, 1, 2, -2, 0, 1, 0, -1, 2, 0, 2, 0, -1, -1, -5, -2, 1, 2, -6, 0, -4, -1],
 [1, 1, -1, 1, 2, -1, 0, 1, 1, 2, 0, -1, -2, 0, -2, 1, 2, 1, -8, 2, 0, 0, -1, -1, -1, -2, -1, 0, -2, -2, -8, 1],
 [2, 2, 2, 2, -2, 2, 0, 2, 2, -2, -6, 2, 0, 0, -2, 2, -8, 2, -2, -2, 0, -6, 2, 2, 2, 0, 2, 0, 0, -2, 4, 2]];
lmf_level := 138;
lmf_primes := [5, 7, 11, 13, 17, 19, 29, 31];

// [1] the transcription: the traces at 2, 3, 23 against the Atkin-Lehner eigenvalues, and
// multiplicativity, which rebuilds all 32 traces of a dimension-1 form from its values at primes
for i in [1..#lmf_label] do
    for j in [1..3] do
        p := [2, 3, 23][j];
        error if lmf_tr[i][p] ne -lmf_dim[i] * lmf_al[i][j],
            Sprintf("LMFDB138: %o has a_%o = %o but Atkin-Lehner eigenvalue %o; a_p = -dim * w_p "
                    * "at a prime exactly dividing the level", lmf_label[i], p, lmf_tr[i][p], lmf_al[i][j]);
    end for;
    if lmf_dim[i] ne 1 then continue; end if;
    lmf_a := [Integers()| 1 : n in [1..32]];
    for n in [2..32] do
        fa := Factorization(n);
        if #fa eq 1 then
            p := fa[1][1]; e := fa[1][2];
            if lmf_level mod p eq 0 then
                lmf_a[n] := lmf_tr[i][p]^e;
            else
                lmf_a[n] := lmf_tr[i][p] * lmf_a[n div p] - (e ge 2 select p * lmf_a[n div p^2] else 0);
            end if;
        else
            lmf_a[n] := &*[lmf_a[fa[j][1]^fa[j][2]] : j in [1..#fa]];
        end if;
    end for;
    error if lmf_a ne lmf_tr[i],
        Sprintf("LMFDB138: the traces of %o are not multiplicative, so the transcription is "
                * "wrong: fetched %o, rebuilt from the primes %o", lmf_label[i], lmf_tr[i], lmf_a);
end for;

// [2] the forms selected by W, under the Jacquet-Langlands sign rule
function lmf_eps(i, w)
    e := 1;
    for j in [1..3] do
        p := [2, 3, 23][j];
        if w mod p eq 0 then e *:= (p eq 23 select lmf_al[i][j] else -lmf_al[i][j]); end if;
    end for;
    return e;
end function;

function lmf_selected(W)
    return [i : i in [1..#lmf_label] | forall{w : w in W | lmf_eps(i, w) eq 1}];
end function;

// [3] the curves: every stored entry of genus at least 1, and the top curve rebuilt
P<x> := PolynomialRing(Rationals());
lmf_models := eval (Read("data/models/models_6_23.m") cat "\nreturn models;");

lmf_tag := []; lmf_W := []; lmf_g := []; lmf_fs := [* *];
for k in Keys(lmf_models) do
    for e in lmf_models[k] do
        if e[1] lt 1 then continue; end if;
        W := {Integers()| a : a in k};
        if Type(e[2]) eq MonStgElt then
            fs := CRVDoubleCovers(e[3]);
        else
            fs := [e[2] + (#e ge 3 select e[3]^2/4 else 0)];     // y^2 + h y = f is y'^2 = f + h^2/4
        end if;
        Append(~lmf_tag, Sprintf("W = %o", Sort(SetToSequence(W))));
        Append(~lmf_W, W); Append(~lmf_g, e[1]); Append(~lmf_fs, fs);
    end for;
end for;

// the top curve as the fibre product of three stored double covers over the star line
lmf_poly := AssociativeArray();
for k in Keys(lmf_models) do
    if #lmf_models[k] gt 0 and Type(lmf_models[k][1][2]) ne MonStgElt then
        lmf_poly[{Integers()| a : a in k}] := lmf_models[k][1][2]
            + (#lmf_models[k][1] ge 3 select lmf_models[k][1][3]^2/4 else 0);
    end if;
end for;
lmf_ups := {U : U in AtkinLehnerDoubleCoversOver({Integers()|1}, 6, 23) | IsDefined(lmf_poly, U)};
lmf_ok, lmf_gens := FibreProductGenerators({Integers()|1}, lmf_ups, 6, 23);
error if not lmf_ok, "LMFDB138: the stored double covers of X_0(6,23)* do not generate the top curve";
Append(~lmf_tag, "W = [ 1 ] rebuilt from the double covers");
Append(~lmf_W, {Integers()|1}); Append(~lmf_g, 5);
Append(~lmf_fs, [lmf_poly[U] : U in lmf_gens]);

// [4] dimensions against genera, and point counts against the newform traces
lmf_bad := []; lmf_n1 := 0; lmf_n2 := 0;
for i in [1..#lmf_tag] do
    sel := lmf_selected(lmf_W[i]);
    error if &+[Integers()| lmf_dim[j] : j in sel] ne lmf_g[i],
        Sprintf("LMFDB138: %o has genus %o but the newforms fixed by W have total dimension %o (%o)",
                lmf_tag[i], lmf_g[i], &+[Integers()| lmf_dim[j] : j in sel], [lmf_label[j] : j in sel]);
    rational := forall{j : j in sel | lmf_dim[j] eq 1};
    good := 0;
    for p in lmf_primes do
        Kp := RationalFunctionField(GF(p)); L := Kp;
        built := true;
        try
            for f in lmf_fs[i] do
                RY<Y> := PolynomialRing(L);
                L := FunctionField(Y^2 - L!Evaluate(PolynomialRing(GF(p))!f, Kp.1));
            end for;
        catch err built := false;
        end try;
        if not built or Genus(L) ne lmf_g[i] then continue; end if;
        good +:= 1;
        n1 := #Places(L, 1);
        e1 := p + 1 - &+[Integers()| lmf_tr[j][p] : j in sel];
        if n1 ne e1 then
            Append(~lmf_bad, Sprintf("%o over F_%o: %o points, the newforms %o give %o",
                                     lmf_tag[i], p, n1, [lmf_label[j] : j in sel], e1));
        end if;
        lmf_n1 +:= 1;
        if not (rational or p^2 le 32) then continue; end if;
        // Tr a_{p^2} = a_p^2 - p for a form of dimension 1; otherwise it is the fetched trace
        tr2 := &+[Integers()| lmf_dim[j] eq 1 select lmf_tr[j][p]^2 - p else lmf_tr[j][p^2] : j in sel];
        n2 := n1 + 2 * #Places(L, 2);
        e2 := p^2 + 1 - (tr2 - p * &+[Integers()| lmf_dim[j] : j in sel]);
        if n2 ne e2 then
            Append(~lmf_bad, Sprintf("%o over F_%o^2: %o points, the newforms give %o",
                                     lmf_tag[i], p, n2, e2));
        end if;
        lmf_n2 +:= 1;
    end for;
    error if good lt 6,
        Sprintf("LMFDB138: %o has good reduction at only %o of the %o primes, too few to check",
                lmf_tag[i], good, #lmf_primes);
end for;

error if #lmf_tag lt 14 or lmf_n1 eq 0,
    Sprintf("LMFDB138: NO EVIDENCE -- %o curves and %o point counts", #lmf_tag, lmf_n1);
error if not IsEmpty(lmf_bad),
    Sprintf("LMFDB138: %o disagreement(s) with the newforms of level 138: %o", #lmf_bad, lmf_bad);

printf "  ok (%o curves, %o counts over F_p and %o over F_p^2)\n", #lmf_tag, lmf_n1, lmf_n2;
