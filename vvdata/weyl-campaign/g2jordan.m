// The 2-part of Dbar = (A, -q) across the two parities (2|N vs 2|D): its Jordan
// type and the Gauss function G2(c) = \mathscr{G}(c, 0; Dbar_2) for odd c --
// closing step (i) of the theta_g derivation (fitted: G2a(x) = (2/x) e_8(x) on
// 2|N bases, G2b = -G2a on 2|D bases).
//
//   magma -b g2jordan.m
AttachSpec("ShimuraQuotients.spec");

PREC := 60;
CC := ComplexField(PREC); ii := CC.1; pi := Pi(CC);
ee := func< z | Exp(2*pi*ii*z) >;
R8 := RealField(8);

for DN in [ <15,2>, <21,2>, <33,2>, <35,2>, <10,7>, <6,11>, <6,5>, <22,3>, <34,5>, <14,3> ] do
    D := DN[1]; N := DN[2];
    Ld := ShimuraCurveLattice(D, N);
    G := Ld`disc_grp;
    Qr := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
    elts := [ g : g in G ];
    n := #elts;
    frac := func< r | r - Floor(r) >;
    Qb := AssociativeArray();
    for g in elts do
        v := ChangeRing(g@@Ld`to_disc, Rationals());
        Qb[g] := frac(-(v*Qr, v)/(2*dn^2));
    end for;
    // the 2-part: elements killed by the 2-power part of the exponent
    ex := LCM([ m : m in Moduli(G) | m gt 1 ]);
    e2 := 2^Valuation(ex, 2);
    odd := ex div e2;
    two := [ g : g in elts | IsZero(e2*g) ];
    qvals := Sort([ Qb[g] : g in two ]);
    // G2(c) for odd c: (1/sqrt(|J| |J[c]|)) sum e(c Qbar) over the 2-part; |J[c]| = |J|
    printf "X0^%o(%o): |2-part| = %o, Q values %o\n", D, N, #two, qvals;
    printf "   parity: %o | %o\n", IsEven(D) select "2|D" else "2|N", D*N;
    for c in [1, 3, 5, 7] do
        S := &+[ CC | ee(frac(c*Qb[g])) : g in two ];
        g2 := S / CC!#two;
        pred := KroneckerSymbol(2, c) * ee(CC!c/8);
        r := g2/pred;
        printf "   G2(%o) arg8 = %o   G2/[ (2/c)e_8(c) ] = (%o, %o)\n", c,
            R8!(8*Arg(g2)/(2*Pi(RealField(20)))), R8!Re(r), R8!Im(r);
    end for;
end for;
quit;
