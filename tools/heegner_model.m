// A model of X = X_0^D(N)/W (genus 3) from the Heegner heights of its branch points on the elliptic
// quotient E = X/<W, w_U>, by ModelByHeegnerHeights, with the result reduced by SmallHyperellipticModel
// and, when an involution has rational fixed points, KleinFourEvenModel.
//
//   magma -b D:=1326 N:=1 W:=1,17,78,1326 U:=6 tools/heegner_model.m < /dev/null
//   magma -b D:=210 N:=13 W:=1,2,21,42,65,130,1365,2730 U:=5 tools/heegner_model.m < /dev/null
//
// Every curve of the isogeny class of E is tried; on all but the quotient itself the predicted heights are
// not heights of lattice points and the search returns nothing. The run prints the surviving
// configurations, the twist, and the model; the intrinsic omits configurations that are rational
// translates or negatives of an earlier one, and a survivor isomorphic to a model already printed is
// marked as such.
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
// For a cover whose canonical conic has no rational point there is no model y^2 = f(x), and the model
// over the conic is printed instead: a X1^2 + b X2^2 + c X3^2 = 0 with y^2 = A(X1, X2). Over the line
// s = X1/X2 such a curve is the compositum of Z^2 = -(a s^2 + b)/c and y^2 = A(s, 1), so an
// isomorphism of two of them descends to a Mobius map of that line, and the map is found from the
// isomorphism of the genus-2 quotients v^2 = -c (a s^2 + b) A(s, 1), where Magma's IsIsomorphic is
// usable. The two conditions below say that the map carries one conic to the other and one quartic to
// a square multiple of the other, which is what makes the two curves isomorphic over Q.
PS<Sv> := PolynomialRing(Rationals());
quartic := func<A, R | PS![MonomialCoefficient(A, R.1^i*R.2^(4-i)) : i in [0..4]]>;
function same_conic_model(m1, m2)
    a1, b1, c1, A1 := Explode(m1); a2, b2, c2, A2 := Explode(m2);
    C1 := HyperellipticCurve(-c1*(a1*Sv^2 + b1)*A1);
    C2 := HyperellipticCurve(-c2*(a2*Sv^2 + b2)*A2);
    if not IsIsomorphic(C1, C2) then return false; end if;
    _, phi := IsIsomorphic(C1, C2);
    eqs := DefiningPolynomials(phi);
    num := Evaluate(eqs[1], [Sv, 0, 1]); den := Evaluate(eqs[3], [Sv, 0, 1]);
    q := a2*num^2 + b2*den^2;
    if Degree(q) ne 2 or q ne Coefficient(q, 2)*(a1*Sv^2 + b1) then return false; end if;
    rho := Coefficient(q, 2);
    A2mu := &+[Coefficient(A2, i)*num^i*den^(4-i) : i in [0..4]];
    if A2mu ne (LeadingCoefficient(A2mu)/LeadingCoefficient(A1))*A1 then return false; end if;
    return IsSquare(rho*c1*c2) and IsSquare(LeadingCoefficient(A2mu)/LeadingCoefficient(A1));
end function;
conicmodels := [];
for E in cls do
    printf "\n== E = %o ==\n", CremonaReference(E);
    surv := ModelByHeegnerHeights(X, Wp, E);
    for s in surv do
        printf "rational parts %o, cover %o twisted by %o, noether %o\n", s[1], s[3], s[4], s[5];
        if s[5] ne 5 then printf "   not hyperelliptic\n"; continue; end if;
        if Type(s[6]) ne CrvHyp then
            c := TwistedCover(DoubleCoversWithBranchDivisor(E, &+[Divisor(pl) : pl in {@ Place(P) : P in s[2] @}])[s[3]], s[4]);
            conic, G := HyperellipticModelOverConic(c);
            Rc := Parent(conic);
            printf "   hyperelliptic, canonical conic without a rational point:\n      %o = 0,  y^2 = %o\n", conic, G;
            quad := [MonomialCoefficient(conic, Rc.i^2) : i in [1..3]];
            cm := <quad[1], quad[2], quad[3], quartic(G, Rc)>;
            if exists{m : m in conicmodels | same_conic_model(m, cm)} then
                printf "      (the same curve as a model above)\n";
            else
                Append(~conicmodels, cm);
            end if;
            continue;
        end if;
        H := SmallHyperellipticModel(s[6]);
        ok, He := KleinFourEvenModel(H);
        if ok and #Sprint(HyperellipticPolynomials(He)) lt #Sprint(HyperellipticPolynomials(H)) then H := He; end if;
        f := HyperellipticPolynomials(H);
        new := not exists{m : m in models | IsIsomorphic(m, H)};
        printf "   model y^2 = %o%o\n", f, new select "" else "   (isomorphic to a model above)";
        if new then Append(~models, H); end if;
    end for;
end for;
printf "\n%o model(s) y^2 = f(x) and %o over a conic without rational points, up to Q-isomorphism\n",
    #models, #conicmodels;
