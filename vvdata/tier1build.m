// Tier 1 of the post-m0 backlog sweep: build the sub-hyperelliptic cover
// models of a previously-failing base on merged main (exact m=0 term live).
//
//   magma -b DD:=15 NN:=2 vvdata/tier1build.m
// Run from the repo root (needs the polymake/ cache directory).
AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 1);

D := 15; N := 2;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;

t0 := Cputime();
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
printf "BASE %o %o: %o covers, star genus %o\n", D, N, #star`CoveredBy, star`g;

all_eqns, all_ws := AllEquationsAboveCovers(star, curves);

printf "==== EQUATIONS X0^%o(%o) ====\n", D, N;
for k in Keys(all_eqns) do
    printf "KEY %o (curve W = %o, genus %o):\n", k, curves[k]`W, curves[k]`g;
    inner := all_eqns[k];
    for k2 in Keys(inner) do
        printf "  [%o] %o\n", k2, inner[k2];
    end for;
end for;
printf "==== AL DATA ====\n";
for k in Keys(all_ws) do
    printf "KEY %o:\n", k;
    inner := all_ws[k];
    for k2 in Keys(inner) do
        printf "  [%o] %o\n", k2, inner[k2];
    end for;
end for;

// models_D_N.m-formatted block (the data/models/README.md conventions)
printf "==== MODELS FILE ====\n";
printf "// Subhyperelliptic cover models for X_0(%o,%o)* -- Guo-Yang / AllEquationsAboveCovers\n", D, N;
printf "// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).\n";
printf "P<x> := PolynomialRing(Rationals());\n";
printf "models := AssociativeArray();\n";
for k in Keys(all_eqns) do
    inner := all_eqns[k];
    entries := [ ];
    for k2 in Keys(inner) do
        crv := inner[k2];
        if Type(crv) eq CrvHyp then
            f, h := HyperellipticPolynomials(crv);
            Append(~entries, Sprintf("<%o, P!%o, P!%o>", curves[k]`g,
                Coefficients(f), Coefficients(h)));
        else
            Append(~entries, Sprintf("<%o, \"CRV\", %o>", curves[k]`g,
                [ Sprintf("%o", p) : p in DefiningPolynomials(crv) ]));
        end if;
    end for;
    if #entries gt 0 then
        printf "models[%o] := [* %o *];\n",
            Sort(SetToSequence(curves[k]`W)), Join(entries, ", ");
    end if;
end for;
printf "==== END MODELS FILE ====\n";
printf "TIER1 %o_%o DONE in %o s\n", D, N, Cputime(t0);
quit;
