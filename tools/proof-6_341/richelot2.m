SetColumns(0);
P<x> := PolynomialRing(Rationals());
C := HyperellipticCurve(P![9,-33,9,29,-12,6], P![0,1,1]);
Cs := SimplifiedModel(C); J := Jacobian(Cs);
seen := [Cs]; queue := [Cs];
while #queue gt 0 do
    X := queue[1]; queue := queue[2..#queue];
    for A in RichelotIsogenousSurfaces(Jacobian(X)) do
        if Type(A) ne JacHyp then printf "  product of elliptic curves: %o\n", A; continue; end if;
        Y := SimplifiedModel(ReducedMinimalWeierstrassModel(Curve(A)));
        if exists{Z : Z in seen | IsIsomorphic(Z, Y)} then continue; end if;
        Append(~seen, Y); Append(~queue, Y);
    end for;
end while;
printf "%o curves in the (2,2)-isogeny graph over Q:\n", #seen;
for i->Y in seen do
    Ym := ReducedMinimalWeierstrassModel(Y); Jy := Jacobian(Y);
    pts := Points(Y : Bound := 5000);
    printf "[%o] %o\n     disc primes %o; conductor %o\n     #points(h<=5000) %o, fibres %o, Weierstrass %o\n     torsion %o, Jac[2] order %o, rank bounds %o\n",
        i, Ym, PrimeDivisors(Integers()!Discriminant(Y)), Factorisation(Conductor(Ym)),
        #pts, Sort(Setseq({ pt[1]/pt[3] : pt in pts | pt[3] ne 0 and pt[2] ne 0 })), [pt : pt in pts | pt[2] eq 0],
        Invariants(TorsionSubgroup(Jy)), #TwoTorsionSubgroup(Jy), [RankBounds(Jy)];
    printf "     counts agree with C at p<60: %o\n", forall{p : p in PrimesUpTo(60) | 2046 mod p ne 0 and #Points(ChangeRing(Y, GF(p))) eq #Points(ChangeRing(Cs, GF(p)))};
end for;
quit;
