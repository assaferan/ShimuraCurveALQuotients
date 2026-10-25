SetColumns(0);
P<x> := PolynomialRing(Rationals());
C := HyperellipticCurve(P![9,-33,9,29,-12,6], P![0,1,1]); Cs := SimplifiedModel(C);
a := [ 446378544780376706967613607/10571649294837312, 28355392880146306419804791/169146388717396992, 34561220587056945266945183/676585554869587968 ];
IC := [ Rationals() | 1, a[2]/a[1], a[3]/a[1], 1/a[1] ];
CA := ReducedMinimalWeierstrassModel(HyperellipticCurveFromIgusaClebsch(IC));
// strip the twist: the discriminant primes outside {2,3,11,31} appear to the power 10 = the weight of a quadratic twist
big := &*[ Integers() | p : p in PrimeDivisors(Integers()!Discriminant(CA)) | p notin {2,3,11,31} ];
printf "twisting by %o\n", big;
CA := ReducedMinimalWeierstrassModel(QuadraticTwist(CA, big));
printf "after: disc %o\n   %o\n", Factorisation(Integers()!Discriminant(CA)), CA;
good := [p : p in PrimesUpTo(80) | 2046 mod p ne 0];
function cnt(X) Xs := SimplifiedModel(X); out := []; for p in good do try n := #Points(ChangeRing(Xs, GF(p))); catch e n := -1; end try; Append(~out, n); end for; return out; end function;
target := cnt(Cs);
for d in [1, -1, 2, -2, 3, -3, 6, -6, 11, -11, 22, -22, 33, -33, 66, -66, 31, -31, 62, -62, 93, -93, 186, -186, 341, -341, 682, -682, 1023, -1023, 2046, -2046] do
    Xd := ReducedMinimalWeierstrassModel(QuadraticTwist(CA, d)); c := cnt(Xd);
    if c eq target then
        printf "TWIST d = %o matches the counts of C at all %o good primes below 80\n", d, #good;
        CA := Xd; break;
    end if;
end for;
printf "C_A : %o\n", CA;
printf "isomorphic to C: %o\n", IsIsomorphic(SimplifiedModel(CA), Cs);
printf "discriminant primes %o, conductor %o\n", PrimeDivisors(Integers()!Discriminant(CA)), Factorisation(Conductor(CA));
CAs := SimplifiedModel(CA); JA := Jacobian(CAs);
printf "simplified: y^2 = %o\n   factors %o\n", HyperellipticPolynomials(CAs), Factorisation(HyperellipticPolynomials(CAs));
printf "Jac(C_A)(Q): 2-torsion order %o, torsion %o, rank bounds %o\n", #TwoTorsionSubgroup(JA), Invariants(TorsionSubgroup(JA)), [RankBounds(JA)];
pts := Points(CAs : Bound := 5000);
printf "%o rational points of height <= 5000; fibres %o; Weierstrass %o\n", #pts, Sort(Setseq({ pt[1]/pt[3] : pt in pts | pt[3] ne 0 and pt[2] ne 0 })), [pt : pt in pts | pt[2] eq 0];
printf "Richelot-isogenous surfaces of Jac(C_A) over Q:\n";
for B in RichelotIsogenousSurfaces(JA) do
    if Type(B) eq JacHyp then Y := ReducedMinimalWeierstrassModel(Curve(B)); printf "   %o   (isomorphic to C: %o, to C_A: %o)\n", Y, IsIsomorphic(SimplifiedModel(Y), Cs), IsIsomorphic(SimplifiedModel(Y), CAs); else printf "   %o\n", B; end if;
end for;
printf "CA_MODEL := <%o, %o>;\n", Coefficients(HyperellipticPolynomials(CA)), Coefficients(HyperellipticPolynomials(CA, 2));
quit;
