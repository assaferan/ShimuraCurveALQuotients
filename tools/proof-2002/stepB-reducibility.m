// STEP B: which primes lambda of O_K could carry a Q-rational isogeny out of A_f?
//
// A Q-isogeny of A_f that is not an endomorphism needs a Galois-stable subgroup of some A_f[lambda]
// not of the form A_f[ideal], hence a REDUCIBLE mod-lambda representation (or an isomorphism
// rho_lambda = rho_lambda' at a split prime, treated separately).  If rho_lambda is reducible its
// semisimplification is chi_1 + chi_2 with chi_1 chi_2 the mod-l cyclotomic character and both
// unramified outside 2002 l.  Characters with values in the residue field have order prime to l, so
// there is no wild l-ramification and the l-part of each is a power of the cyclotomic character;
// and the conductors of chi_1, chi_2 multiply to a divisor of 2002, so each has conductor dividing
// 2002 -- hence dividing 1001, there being no nontrivial character of conductor 2.  At a prime
// p = 1 mod 2002 both tame parts are trivial, and weight two leaves a_p = 1 + p (mod lambda).
// So  lambda | gcd over such p of Norm(a_p - 1 - p),  with NO assumption on the characters: an
// exact, finite bound, after which the survivors are checked one by one.
SetColumns(0);
MS := ModularSymbols(2002, 2); NS := NewSubspace(CuspidalSubspace(MS));
F := 0;
for G in NewformDecomposition(NS) do
    if Dimension(G) ne 4 then continue; end if;
    if [ Integers() ! (Eltseq(AtkinLehnerOperator(G, p))[1]) : p in [2,7,11,13] ] eq [-1,-1,-1,-1]
        then F := G; break; end if;
end for;
PQ<x> := PolynomialRing(Rationals());
// a_p lives in K = Q(sqrt 5); Norm(a_p - c) = Resultant of (x^2 - t x + n) with (x - c)
function NormOfApMinus(p, c)
    cp := CharacteristicPolynomial(HeckeOperator(F, p));
    ok, q := IsSquare(cp); assert ok;         // cp = (x^2 - t x + n)^2 always
    return Integers() ! Evaluate(q, c);       // Norm(a_p - c) = prod over conjugates
end function;
ps := [ p : p in PrimesInInterval(2003, 200000) | p mod 2002 eq 1 ];
printf "%o primes p = 1 mod 2002 below 200000\n", #ps;
g := 0;
printf "\n%-8o %-26o %o\n", "p", "Norm(a_p - 1 - p)", "running gcd";
for p in ps do
    v := NormOfApMinus(p, 1 + p);
    g := GCD(g, v);
    printf "%-8o %-26o %o\n", p, v, g;
end for;
printf "\ngcd over all of them: %o = %o\n", g, Factorisation(g);
printf "so a reducible rho_lambda needs lambda | %o -- these, and the primes above 2, 7, 11, 13,\n", g;
printf "are the only candidates for a nontrivial Q-isogeny out of A_f, and each is now a finite check.\n";
quit;
