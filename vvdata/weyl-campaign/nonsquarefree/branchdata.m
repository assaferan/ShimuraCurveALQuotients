// Branch data of the Galois S_3-cover X_0^D(4N')/W_D -> X_0^D(N')^*, N' odd, from counts at the
// base level only.  Over a star point P the fibre pattern is given by the image H_P in S_3 of
// the stabiliser of a point above P: an order-3 elliptic point contributes a 3-cycle, an order-2
// one a transposition, and a w_q-fixed CM point a transposition unless its order is the
// maximal order of Q(sqrt(-q)) (then sqrt(-q) = 1 mod 2 and the image is trivial).  Each star
// point contributes 6 - 6/|H_P| to the ramification, and Riemann-Hurwitz for the degree-6 map
// between a genus-0 base and X/W_D demands R = 10 + 2 g(X/W_D).
AttachSpec("ShimuraQuotients.spec");
hall := func<M | [d : d in Divisors(M) | GCD(d, M div d) eq 1]>;
for b in [[15,1],[21,1],[33,1]] do
    D := b[1]; Np := b[2]; N := 4*Np;
    W := hall(D*Np);
    gWD := GenusShimuraCurveQuotient(D, N, Set(hall(D)));
    gstar := GenusShimuraCurveQuotient(D, Np, Set(W));
    printf "\n== X_0^%o(%o)/W_%o (genus %o) -> X_0^%o(%o)^* (genus %o): need R = %o\n", D, N, D, gWD, D, Np, gstar, 10 + 2*gWD;
    // discriminants of the special points: elliptic, and fixed by some w_q
    discs := {-3, -4};
    fixes := AssociativeArray();          // q -> (disc -> fixed count)
    for q in W do
        if q eq 1 then continue; end if;
        fixes[q] := NumFixedPointsByCMOrder(D, Np, q);
        discs join:= Set(Keys(fixes[q]));
    end for;
    R := 0;
    for d in Sort(SetToSequence(discs)) do
        Rd := QuadraticOrder(BinaryQuadraticForms(d));
        n := NumberOfOptimalEmbeddings(Rd, D, Np);
        if n eq 0 then continue; end if;
        fixv := [n] cat [IsDefined(fixes[q], d) select fixes[q][d] else 0 : q in W | q ne 1];
        orbits := &+fixv / #W;
        stab := #W * orbits / n;
        fixers := [q : q in W | q ne 1 and IsDefined(fixes[q], d) and fixes[q][d] eq n];
        // image in S_3 of each fixing w_q: trivial iff d = -q with q = 3 mod 4 (maximal order)
        transp := [q : q in fixers | not (d eq -q)];
        e := (d eq -3) select 3 else ((d eq -4) select 2 else 1);
        if e eq 3 then
            H := (#transp gt 0) select 6 else 3;
        elif e eq 2 then
            H := 2; note := (#transp gt 0) select " (w-transpositions too: H could be S_3)" else "";
        else
            H := (#transp gt 0) select 2 else 1; note := (#transp gt 1) select " (two w-transpositions: H could be S_3)" else "";
        end if;
        if e eq 3 then note := ""; end if;
        contrib := orbits * (6 - 6/H);
        R +:= contrib;
        printf "  d = %-5o points %-3o star points %-3o stab %-2o e=%o fixers %-12o |H| = %o  contributes %o%o\n",
               d, n, orbits, stab, e, fixers, H, contrib, note;
    end for;
    ok := (R eq 10 + 2*gWD) select "Riemann-Hurwitz OK" else "MISMATCH";
    printf "  total R = %o   (%o)\n", R, ok;
end for;
