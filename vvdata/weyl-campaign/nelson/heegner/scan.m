// Calibration points for the height formula: every stored genus-3 hyperelliptic quotient
// y^2 + h y = f whose completed octic F = 4f + h^2 is even (involution x -> -x) or palindromic
// (x -> 1/x), with the elliptic quotient E having a rational point. For each branch point: the
// field of definition, the canonical height, and the ratio to L(E,1) L'(E^d,1) sqrt|d| / Area.
AttachSpec("ShimuraQuotients.spec"); SetColumns(0); SetDefaultRealFieldPrecision(30);
P<x> := PolynomialRing(Rationals());
files := Split(Pipe("ls data/models", ""), "\n");
nfound := 0;
for fn in files do
    ok, _, m := Regexp("models_([0-9]+)_([0-9]+).m", fn);
    if not ok then continue; end if;
    D := StringToInteger(m[1]); N := StringToInteger(m[2]);
    models := eval (Read("data/models/" cat fn) cat " return models;");
    for W in Keys(models) do
        for ent in models[W] do
            if ent[1] ne 3 or Type(ent[2]) ne RngUPolElt then continue; end if;
            F := 4*ent[2] + ent[3]^2;
            if Degree(F) ne 8 then continue; end if;
            cases := [];
            if &and[Coefficient(F, i) eq 0 : i in [1,3,5,7]] then
                q := P![Coefficient(F, 2*i) : i in [0..4]];
                Append(~cases, <"even", q, [* <0, 1>, <1, 0> *]>);   // branch fibres over u = 0, oo  (u = x^2)
            end if;
            for cs in cases do
                if cs[1] ne "even" then continue; end if;   // palindromic handled below once it is written
                q := cs[2];
                den := LCM([Denominator(c) : c in Coefficients(q)]);
                C := HyperellipticCurve(q*den^2);
                pts := Points(C : Bound := 3000);
                if #pts eq 0 then printf "%o_%o W=%o: E has no small rational point\n", D, N, W; continue; end if;
                E0, mp := EllipticCurve(C, pts[1]); E, mn := MinimalModel(E0);
                u := Pullback(Inverse(mp), FunctionField(C).1);
                Z, Pl := SignDecomposition(Divisor(u));
                per := Periods(E); A := Abs(Imaginary(ComplexConjugate(per[1])*per[2]));
                L1 := Evaluate(LSeries(E), 1);
                nfound +:= 1;
                printf "%o_%o W=%o: E = %o rank %o\n", D, N, W, CremonaReference(E), Rank(E);
                for pl in Support(Z + Pl) do
                    R := RepresentativePoint(pl); K := Ring(Parent(R));
                    if Type(K) eq FldRat then
                        printf "   rational branch point, height %o\n", Height(mn(E0!R)); continue;
                    end if;
                    d := Discriminant(Integers(AbsoluteField(K)));
                    PK := BaseChange(E, K)!Eltseq(mn(E0(K)!R));
                    h := Height(PK);
                    if d gt 0 then printf "   d = %o (real): height %o\n", d, h; continue; end if;
                    tw := QuadraticTwist(E, d);
                    Ld1 := Evaluate(LSeries(tw), 1 : Derivative := 1);
                    base := L1 * Ld1 * Sqrt(Abs(d)) / A;
                    printf "   d = %o: height %o, L(E,1) = %o, L'(E^d,1) = %o, ratio = %o\n", d, h, L1, Ld1, base eq 0 select "inf" else h/base;
                end for;
            end for;
        end for;
    end for;
end for;
printf "cases with E: %o\nDONE\n", nfound;
