// WHICH curve is the target?  A map out of a quotient X_0^D(N)/U down to X^*(D,N') exists only if U
// lies in <Gamma_0(N'), W'>, W' = Hall(D N').  The new Atkin-Lehner involutions of the top level
// (w_4, w_9, w_8, ...) do NOT: w_4, say, has reduced norm 4, so it would have to be 2 times a
// norm-one unit of the maximal order, i.e. w_4/2 integral, which a primitive element of norm 4 is
// not -- and classically w_4 = (0 -1; 4 0)/2 is visibly not in SL_2(Z).  So there is NO map
// X^*(D,N) -> X^*(D,N'); the largest quotient that maps down is
//     Y := X_0^D(N) / W',    W' = Hall(D N'),   of degree psi(N)/psi(N') over X^*(D,N').
// Tu's t_4 lives on exactly this Y for 15_4 (he quotients by W_15, not by Hall(60)).  So the
// question that decides whether a HAUPTMODUL exists at all is the genus of Y.
AttachSpec("ShimuraQuotients.spec");
hall := func<M | {d : d in Divisors(M) | GCD(d, M div d) eq 1}>;
psi := func<M | M * &*[Rationals() | 1 + 1/q : q in PrimeDivisors(M)]>;
printf "%-8o %-6o %-26o %-5o %-7o %-9o %o\n", "base", "N'", "U = Hall(DN) cap Hall(DN')", "deg", "g(X/U)", "g(base/U)", "note";
for b in [[15,4,1],[21,4,1],[33,4,1],[15,8,2],[10,9,1],[14,9,1],[22,9,1]] do
    D := b[1]; N := b[2]; Np := b[3];
    Wp := hall(D*Np) meet hall(D*N);   // only these act on BOTH curves
    gY := GenusShimuraCurveQuotient(D, N, Wp);
    gb := GenusShimuraCurveQuotient(D, Np, Wp);
    deg := Integers() ! (psi(N)/psi(Np));
    printf "%-8o %-6o %-26o %-5o %-7o %-9o %o\n", Sprintf("%o_%o", D, N), Np,
           Sprintf("%o", Sort(SetToSequence(Wp))), deg, gY, gb,
           gY eq 0 select "genus 0: a hauptmodul exists" else "higher genus: no hauptmodul";
end for;
