// Compute Kudla's Diff(m) PROPERLY: the places where the ternary quadratic SPACE fails to represent m.
//
// This replaces the crude proxy "#primes where the code's local density vanishes", which was REFUTED on
// X0^15(2). That proxy conflated two different things: the code's W_p is a LATTICE density (it can
// vanish for congruence reasons at p | level even when the space represents m perfectly well), whereas
// Diff is about the SPACE over Q_p, and it must also include the archimedean place.
//
// Criterion (Serre, A Course in Arithmetic, Ch IV Thm 6, rank 3): a nondegenerate ternary form over Q_p
// with discriminant d and Hasse invariant eps represents a in Q_p^* UNLESS
//        a = -d in Q_p^*/(Q_p^*)^2   AND   eps != (-1,-d)_p .
// Hasse invariant for f ~ <a1,a2,a3> is eps_p = prod_{i<j} (a_i,a_j)_p.
//
// usage: magma -b DD:=<D> NN:=<N> MMAX:=<bound> diffset.m
// Emits one [DIFF] line per m: the finite places where the space fails, for BOTH sign conventions of
// the Gram (the code's Q, and -Q which memory argues is the Borcherds lattice of signature (2,1)).

AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 0);

D := StringToInteger(DD); N := StringToInteger(NN);
mmax := StringToInteger(MMAX);
Ld := ShimuraCurveLattice(D, N);
Q := ChangeRing(Ld`Q, Rationals());

// diagonalise the Gram over Q
function diag(G)
    n := Nrows(G);
    M := G;
    T := IdentityMatrix(Rationals(), n);
    // simple symmetric Gram-Schmidt (the forms here are nondegenerate)
    for i in [1..n] do
        if M[i][i] eq 0 then
            // find j > i with M[j][j] != 0, or a nonzero off-diagonal to rotate in
            found := false;
            for j in [i+1..n] do
                if M[j][j] ne 0 then
                    SwapRows(~M, i, j); SwapColumns(~M, i, j); found := true; break;
                end if;
            end for;
            if not found then
                for j in [i+1..n] do
                    if M[i][j] ne 0 then
                        AddRow(~M, 1, j, i); AddColumn(~M, 1, j, i); found := true; break;
                    end if;
                end for;
            end if;
            error if not found, "cannot diagonalise";
        end if;
        for j in [i+1..n] do
            c := -M[i][j]/M[i][i];
            AddRow(~M, c, i, j); AddColumn(~M, c, i, j);
        end for;
    end for;
    return [M[i][i] : i in [1..n]];
end function;

// does the ternary space with diagonal `a` represent t over Q_p ?
// eps and the comparison value depend only on (a,p), so they are precomputed by the caller.
function represents(d, eps, cmp, Kp, t)
    // a = -d as square classes  <=>  -t*d is a square in Q_p
    if not IsSquare(Kp ! (-t*d)) then return true; end if;
    return eps eq cmp;
end function;

for sgn in [1, -1] do
    a := diag(sgn*Q);
    dd := a[1]*a[2]*a[3];
    sig := #[x : x in a | x gt 0];
    printf "[SIG] %o %o sgn=%-3o diag=%o det=%o signature=(%o,%o)\n", D, N, sgn, a, dd, sig, 3-sig;
    ramified := PrimeDivisors(D);
    cand := Sort(SetToSequence(Set(PrimeDivisors(2*Numerator(dd)*Denominator(dd))) join Set(ramified)));
    Kps := AssociativeArray(); epss := AssociativeArray(); cmps := AssociativeArray();
    for p in cand do
        Kps[p] := pAdicField(p, 40);
        epss[p] := HilbertSymbol(Rationals()!a[1], Rationals()!a[2], p)
                 * HilbertSymbol(Rationals()!a[1], Rationals()!a[3], p)
                 * HilbertSymbol(Rationals()!a[2], Rationals()!a[3], p);
        cmps[p] := HilbertSymbol(Rationals()!(-1), Rationals()!(-dd), p);
    end for;
    for m in [1..mmax] do
        bad := [p : p in cand | not represents(dd, epss[p], cmps[p], Kps[p], m)];
        // archimedean: a definite form fails to represent m of the wrong sign
        arch := (sig eq 0) or (sig eq 3 and false);   // negative definite fails for m > 0
        if not IsEmpty(bad) or arch then
            printf "[DIFF] %o %o sgn=%-3o m=%-4o finite_fail=%o arch_fail=%o\n", D, N, sgn, m, bad, arch;
        end if;
    end for;
end for;
exit;
