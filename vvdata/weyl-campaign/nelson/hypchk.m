// Re-derive the hyperellipticity of stored CRV models from their function fields by Noether's
// count (the computation of PR #108's IsHyperellipticFunctionField, inlined), so that a control's
// verdict is reproduced here rather than read from a log.
//   magma -b "ITEMS:=46_3:1,3;46_3:1,6;14_11:1,22" hypchk.m < /dev/null
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
MODELDIR := "/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/oldwork/data/models/";

function CRVCovers(eqns)
    R<s, x, y, y1, y2, y3, y4, z> := PolynomialRing(Rationals(), 8);
    es := [eval ("return " cat st cat ";") : st in eqns];
    nocc := [#[1 : E in es | Degree(E, R.j) gt 0] : j in [1..8]];
    fib := [];
    for i in [1..#es] do
        cand := [j : j in [1..8] | nocc[j] eq 1 and Degree(es[i], R.j) eq 2
                 and Degree(es[i] - MonomialCoefficient(es[i], R.j^2)*R.j^2, R.j) le 0];
        assert #cand eq 1;
        Append(~fib, cand[1]);
    end for;
    base := [j : j in [1..8] | nocc[j] gt 0 and j notin fib];
    assert #base eq 2;
    Pt<t> := PolynomialRing(Rationals());
    vals := [Pt!0 : j in [1..8]];
    vals[base[1]] := t; vals[base[2]] := Pt!1;
    fs := [];
    for i in [1..#es] do
        c := MonomialCoefficient(es[i], R.fib[i]^2);
        G := es[i] - c*R.fib[i]^2;
        Append(~fs, -Evaluate(G, vals) / c);
    end for;
    return fs;
end function;

for it in Split(ITEMS, ";") do
    pp := Split(it, ":");
    key := Sort([ StringToInteger(x) : x in Split(pp[2], ",") ]);
    P<x> := PolynomialRing(Rationals());
    models := eval (Read(MODELDIR cat "models_" cat pp[1] cat ".m") cat "\nreturn models;");
    for e in models[key] do
        if Type(e[2]) ne MonStgElt then
            printf "%o %o: stored as y^2 + h y = f, genus %o, degree %o -> hyperelliptic by its model\n", pp[1], key, e[1], Degree(e[2]);
            continue;
        end if;
        t0 := Cputime();
        fs := CRVCovers(e[3]);
        K := FibreProductFunctionField(fs);
        g := Genus(K);
        V, mp := SpaceOfHolomorphicDifferentials(K);
        ws := [mp(V.i) : i in [1..g]];
        phis := [ws[i] / ws[1] : i in [1..g]];
        prods := [phis[i]*phis[j] : j in [i..g], i in [1..g]];
        nrel := Dimension(Relations(prods, Rationals()));
        printf "%o %o: CRV, genus %o (stored %o), quadric relations %o -> %o  (%.1os)\n", pp[1], key, g, e[1], nrel,
            nrel eq (g-1)*(g-2) div 2 select "HYPERELLIPTIC" else (nrel eq (g-2)*(g-3) div 2 select "NOT hyperelliptic" else "??"),
            Cputime(t0);
    end for;
end for;
printf "HYPCHK DONE\n";
quit;
