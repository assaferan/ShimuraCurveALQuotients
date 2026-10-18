// A model of X = X_0^D(N)/W (genus 3) from the Heegner heights of its branch points on the elliptic
// quotient E = X/<W, w_U>, by ModelByHeegnerHeights, with the result reduced by SmallHyperellipticModel
// and, when an involution has rational fixed points, KleinFourEvenModel.
//
//   magma -b D:=1326 N:=1 W:=1,17,78,1326 U:=6 tools/heegner_model.m < /dev/null
//   magma -b D:=210 N:=13 W:=1,2,21,42,65,130,1365,2730 U:=5 tools/heegner_model.m < /dev/null
//
// Every curve of the isogeny class of E is tried; on all but the quotient itself the predicted heights are
// not heights of lattice points and the search returns nothing. The run prints the surviving
// configurations, the twist, and the model; survivors that are rational translates of one another are
// merged by the intrinsic, those differing by other symmetries (negation, conjugation) are isomorphic
// curves and are reported as such.
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
SetVerbose("HeegnerHeights", 1);
D := StringToInteger(D); N := StringToInteger(N); U := StringToInteger(U);
W := {Integers() | StringToInteger(s) : s in Split(W, ",")};
X := CreateShimuraQuot(D, N, W); X`g := GenusShimuraCurveQuotient(D, N, W);
quots := EllipticQuotientsBySigns(D, N, W);
sel := [q : q in quots | U in q[1]];
error if #sel eq 0, "no elliptic quotient through w_U with a newform of level DN";
Wp, Erep := Explode(sel[1]);
DB := CremonaDatabase();
cls := [EllipticCurve(DB, D*N, i, j) : j in [1..NumberOfCurves(DB, D*N, i)]]
       where i := [i : i in [1..NumberOfIsogenyClasses(DB, D*N)] | IsIsogenous(EllipticCurve(DB, D*N, i, 1), Erep)][1];
printf "X_0^%o(%o)/%o, elliptic quotient by w_%o, isogeny class %o\n", D, N, W, U, [CremonaReference(E) : E in cls];
models := [];
for E in cls do
    printf "\n== E = %o ==\n", CremonaReference(E);
    surv := ModelByHeegnerHeights(X, Wp, E);
    for s in surv do
        printf "class b = %o, cover %o twisted by %o, noether %o\n", s[1], s[4], s[5], s[6];
        if s[6] ne 5 then printf "   not hyperelliptic\n"; continue; end if;
        if Type(s[7]) ne CrvHyp then printf "   hyperelliptic, canonical conic without rational point\n"; continue; end if;
        H := SmallHyperellipticModel(s[7]);
        ok, He := KleinFourEvenModel(H);
        if ok and #Sprint(HyperellipticPolynomials(He)) lt #Sprint(HyperellipticPolynomials(H)) then H := He; end if;
        f := HyperellipticPolynomials(H);
        new := not exists{m : m in models | IsIsomorphic(m, H)};
        printf "   model y^2 = %o%o\n", f, new select "" else "   (isomorphic to a model above)";
        if new then Append(~models, H); end if;
    end for;
end for;
printf "\n%o model(s) up to Q-isomorphism\n", #models;
