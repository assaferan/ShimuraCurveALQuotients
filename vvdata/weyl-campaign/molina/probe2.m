// Probe 2.  Measure, on Molina's ten published equations, the data his Steps 2-4 must produce:
//   leading coefficients a_{R_i}; disc(p_{R_i}); the IndexFormEquation index I_i; resultants;
//   and the valuation of disc(R) at 2 and at each p | D.
// Compare the p | D valuations with the sum of thicknesses over the singular points of the
// quotient special fibre (Lemma 4.1 + Sec 10.1), computed from DualGraphData.

AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
P<x> := PolynomialRing(Rationals());

molina := [* <39, {Integers()|1},     3, [-(7*x^4+79*x^3+311*x^2+497*x+277), x^4+9*x^3+29*x^2+39*x+19]>,
             <55, {Integers()|1},     3, [-(3*x^4-32*x^3+130*x^2-237*x+163), x^4-8*x^3+34*x^2-83*x+81]>,
             <35, {Integers()|1,5},   2, [-x, 9*x+4, 4*x+1, 172*x^3+176*x^2+60*x+7]>,
             <51, {Integers()|1,17},  2, [-x, 7*x^3+52*x^2+116*x+68, x-1, x+3]>,
             <57, {Integers()|1,3},   2, [-(x-9), x^3-19*x^2+119*x-249, 7*x^2-104*x+388]>,
             <65, {Integers()|1,13},  2, [-(x^2-3*x+1), 7*x^4-3*x^3-32*x^2+25*x-5]>,
             <65, {Integers()|1,5},   2, [-(x^2+7*x+9), 7*x^4+81*x^3+319*x^2+508*x+268]>,
             <69, {Integers()|1,23},  2, [-x, x+4, 4*x^4-16*x^3+11*x^2+10*x+3]>,
             <85, {Integers()|1,5},   2, [-(3*x^2-41*x+133), x^4-23*x^3+183*x^2-556*x+412]>,
             <85, {Integers()|1,85},  2, [x^2-3*x+1, x^4+x^3-15*x^2+20*x-8]> *];

// Sum of thicknesses over the singular points of the special fibre of X_0^D/W at p | D.
// Edges of the dual graph are the singular points of X_p; thickness = #(O^*/{+-1}) (Lemma 4.1).
// Section 10.1: w_m with p | m swaps the two branches at a fixed edge, so the image is smooth;
// with p not dividing m it fixes both branches and the thickness doubles; edges moved by w_m
// pair up and keep their thickness.
function QuotientThicknesses(D, p, W)
    data := DualGraphData(D, 1, p);
    nE := data[2]; alE := data[4]; uE := data[7];
    thick := [uE[e] div 2 : e in [1..nE]];
    act := AssociativeArray();
    act[1] := [1..nE];
    for t in alE do act[t[1]] := t[2]; end for;
    for m in Sort(Setseq(W)) do
        if IsDefined(act, m) then continue; end if;
        img := [1..nE];
        for l in PrimeDivisors(m) do
            q := l^Valuation(D, l);
            img := [act[q][img[e]] : e in [1..nE]];
        end for;
        act[m] := img;
    end for;
    out := []; seen := {};
    for e in [1..nE] do
        if e in seen then continue; end if;
        seen join:= {act[m][e] : m in W};
        stab := {m : m in W | act[m][e] eq e};
        if exists{m : m in stab | m mod p eq 0} then continue; end if;
        Append(~out, thick[e] * #stab);
    end for;
    return out, Sort([uE[e] div 2 : e in [1..nE]]);
end function;

for e in molina do
    D := e[1]; W := e[2]; g := e[3]; fac := e[4];
    f := &*fac;
    printf "\n=== D = %o, W = %o, genus %o\n", D, Sort(Setseq(W)), g;
    dR := Integers()!Discriminant(f);
    printf "  a_R = %o,  disc(R) = %o,  nu_2 = %o (4(g+1) = %o)\n",
        LeadingCoefficient(f), Factorization(dR), Valuation(dR,2), 4*(g+1);
    for i in [1..#fac] do
        pi := fac[i]; n := Degree(pi); a := LeadingCoefficient(pi);
        if n eq 1 then
            printf "   factor %o: deg 1, a = %o\n", i, a; continue;
        end if;
        dpi := Integers()!Discriminant(pi);
        dM := Discriminant(MaximalOrder(NumberField(pi)));
        I2 := a^(2*n-2)*dpi/dM;
        ok, I := IsSquare(Numerator(I2)*Denominator(I2));
        printf "   factor %o: deg %o, a = %o, disc(p) = %o, disc(M) = %o, index = %o\n",
            i, n, a, Factorization(dpi), dM, ok select I else "NOT A SQUARE";
    end for;
    for i in [1..#fac] do for j in [i+1..#fac] do
        r := Integers()!Resultant(fac[i], fac[j]);
        printf "   Res(%o,%o) = %o = %o\n", i, j, r, Factorization(r);
    end for; end for;
    for p in PrimeDivisors(D) do
        th, xth := QuotientThicknesses(D, p, W);
        printf "  p = %o: nu_p(disc R) = %o ; X_p thicknesses %o ; quotient %o, sum %o\n",
            p, Valuation(dR, p), xth, Sort(th), &+th;
    end for;
end for;

printf "\nPROBE2 DONE\n";
quit;
