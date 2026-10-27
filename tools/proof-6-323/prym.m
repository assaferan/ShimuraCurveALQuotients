// L-polynomials of the Prym of Y = X/H -> E = X/W at small primes: if irreducible of degree 4 over Q at one prime,
// the Prym is simple over Q, so Y is not the twist of the hyperelliptic cover (whose Prym is a product of two
// elliptic curves over Q).
AttachSpec("ShimuraQuotients.spec");
R<x> := PolynomialRing(Rationals());
procedure PrymPolys(D, N, gens, Elab)
    W := {m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1};
    mul := func< a, b | (a*b) div GCD(a, b)^2 >;
    H := {1} join gens; repeat old := H; H := {mul(a, b) : a in H, b in H}; until H eq old;
    Y := CreateShimuraQuot(D, N, H); Y`g := GenusShimuraCurveQuotient(D, N, H); Y`CurveID := 0;
    E := EllipticCurve(Elab);
    printf "%o_%o, H = %o (order %o), genus %o, E = %o\n", D, N, Sort(Setseq(H)), #H, Y`g, Elab;
    for p in [q : q in PrimesInInterval(7, 40) | (D*N) mod q ne 0] do
        // L-polynomial of Jac(Y) from the counts over F_p, F_p^2, F_p^3 (genus 3)
        s := [ p^i + 1 - ComputePointsViaTrace(Y, p, i) : i in [1..3] ];
        e := [Rationals() | 1];
        for k in [1..3] do Append(~e, (&+[ (-1)^(i-1) * e[k-i+1] * s[i] : i in [1..k] ]) / k); end for;
        Ek := func< k | k le 3 select e[k+1] else p^(k-3) * e[6-k+1] >;
        LY := &+[ (-1)^k * Ek(k) * x^k : k in [0..6] ];
        LE := 1 - TraceOfFrobenius(E, p) * x + p * x^2;
        ok, Lprym := IsDivisibleBy(LY, LE);
        assert ok;
        fac := Factorization(Lprym);
        printf "  p = %o: Prym L-polynomial %o, factor degrees %o\n", p, Lprym, [Degree(f[1]) : f in fac];
    end for;
end procedure;
PrymPolys(6, 323, {6, 19, 34}, "1938h1");
quit;
