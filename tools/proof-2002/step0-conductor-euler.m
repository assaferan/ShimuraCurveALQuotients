// STEP 0 OF A PROOF THAT X_0^2002(1)* IS THE CURVE C: two exact necessary conditions.
//   (i)  the conductor of Jac(C) is 2002^2, as for a GL_2-type surface of level 2002;
//   (ii) the FULL Euler factor of Jac(C) at each good prime equals that of A_f, i.e. the product
//        over the two conjugates of (1 - a_p^s T + p T^2).  The earlier check compared only the
//        trace; the degree-four factor also pins Norm(a_p), which the trace does not.
SetColumns(0);
PQ<x> := PolynomialRing(Rationals()); PT<T> := PolynomialRing(Rationals());
C := HyperellipticCurve(16*x^6+4*x^5+180*x^4-456*x^3+196*x^2+64*x+4, x^2+x);
printf "C : %o\n", C;
printf "discriminant: %o\n", Factorisation(Integers() ! Discriminant(C));
cond := Conductor(C);
printf "conductor of Jac(C): %o = %o   (2002^2 = %o)\n", cond, Factorisation(cond), 2002^2;
// the newform side: on the four-dimensional modular-symbol space T_p has characteristic polynomial
// (x^2 - t x + n)^2 with t = Tr a_p and n = Norm a_p, so L_p(T) = 1 - tT + (n+2p)T^2 - tpT^3 + p^2T^4
MS := ModularSymbols(2002, 2); NS := NewSubspace(CuspidalSubspace(MS));
F := 0;
for G in NewformDecomposition(NS) do
    if Dimension(G) ne 4 then continue; end if;
    sg := [ Integers() ! (Eltseq(AtkinLehnerOperator(G, p))[1]) : p in [2,7,11,13] ];
    if sg eq [-1,-1,-1,-1] then F := G; break; end if;    // classical sign -1 at every p | D
end for;
printf "\n%-5o %-34o %-34o %o\n", "p", "L_p of Jac(C)", "L_p of A_f", "equal";
allok := true;
for p in PrimesInInterval(3, 120) do
    if 2002 mod p eq 0 then continue; end if;
    Lc := LPolynomial(ChangeRing(C, GF(p)));
    cp := CharacteristicPolynomial(HeckeOperator(F, p));   // (x^2 - t x + n)^2
    fac := Factorisation(cp);
    q := fac[1][1];                                        // x^2 - t x + n
    t := -Coefficient(q, 1); n := Coefficient(q, 0);
    Lf := 1 - t*T + (n + 2*p)*T^2 - t*p*T^3 + p^2*T^4;
    ok := (PT ! Lc) eq Lf;
    if not ok then allok := false; end if;
    if p le 31 or not ok then printf "%-5o %-34o %-34o %o\n", p, Lc, Lf, ok; end if;
end for;
printf "...\nfull Euler factors agree at every good prime below 120: %o\n", allok;
quit;
