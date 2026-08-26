// fastcosets: a drop-in replacement for VVCosetReps that is O(psi(M)*phi(M))
// instead of O(M^3).
//
// VVCosetReps loops over ALL M^2 pairs (c,d) and, for each one that is
// primitive mod M, builds the whole phi(M)-element orbit under (Z/M)^* just to
// test whether it has been seen.  That is ~M^3 work: measured 0.37 s at
// M = 132, 9.8 s at M = 420, 40.2 s at M = 660 -- clean cubic scaling, so
// M = 4620 would be about 3.8 HOURS before the driver computes anything.
//
// What it actually emits, though, is easy to characterise:
//   * the (c,d) it builds a matrix from is the LEXICOGRAPHIC MINIMUM of its
//     orbit -- when an orbit is first met, the current (c,d) is lex-smallest
//     in it, since any lex-smaller member would already have been marked seen;
//   * the emission ORDER is lex order of those minima, since c is the outer
//     loop and d the inner one.
// So the output is exactly: the lex-min of each orbit, sorted lexicographically,
// pushed through the same lifting/XGCD code.
//
// This version enumerates one representative per orbit directly, by CRT over
// the prime powers: P^1(Z/M) = prod_p P^1(Z/p^e), with the standard local
// representatives (x : 1) for x in Z/p^e and (1 : p*y) for y in Z/p^{e-1}.
// That is psi(M) = prod (p^e + p^{e-1}) points, one per orbit by construction
// -- 13824 at M = 4620 -- and only then is each orbit's lex-min taken.
//
// Both the minimum and the sort are written out explicitly rather than left to
// Magma's tuple comparison, so the result does not depend on how `Min`/`Sort`
// happen to order tuples.  Validated against VVCosetReps for equality as a
// SEQUENCE OF MATRICES (see fastcosets_check.m).

fastCosetReps := function(M)
    // one representative per orbit, via CRT over the prime powers
    cands := [ <0, 0> ];
    modsofar := 1;
    for pe in Factorization(M) do
        p := pe[1]; e := pe[2]; q := p^e;
        L := [ <x, 1> : x in [0..q-1] ] cat [ <1, p*y> : y in [0..p^(e-1) - 1] ];
        nxt := [ ];
        for cd in cands do
            for uv in L do
                Append(~nxt, < CRT([cd[1], uv[1]], [modsofar, q]),
                              CRT([cd[2], uv[2]], [modsofar, q]) >);
            end for;
        end for;
        cands := nxt;
        modsofar *:= q;
    end for;

    units := [ u : u in [1..M] | GCD(u, M) eq 1 ];
    mins := [ ];
    for cd in cands do
        bc := M + 1; bd := M + 1;
        for u in units do
            c2 := (u*cd[1]) mod M; d2 := (u*cd[2]) mod M;
            if c2 lt bc or (c2 eq bc and d2 lt bd) then bc := c2; bd := d2; end if;
        end for;
        Append(~mins, <bc, bd>);
    end for;
    // one candidate per orbit by construction, so the minima are already
    // distinct -- assert it rather than assume it
    error if #Seqset(mins) ne #mins, "CRT enumeration produced a repeated orbit";
    Sort(~mins, func< x, y | x[1] ne y[1] select x[1] - y[1] else x[2] - y[2] >);

    reps := [ ];
    for key in mins do
        cc := key[1]; dd := key[2];
        if cc eq 0 and dd eq 0 then cc := M; dd := 1; end if;
        while GCD(cc, dd) ne 1 do dd +:= M; end while;
        _, a, b := XGCD(dd, -cc);                      // a*dd - b*cc = 1
        Append(~reps, Matrix(Integers(), 2, 2, [a, b, cc, dd]));
    end for;
    return reps;
end function;
