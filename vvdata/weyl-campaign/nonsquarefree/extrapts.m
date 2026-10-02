// The branch points the enumeration missed.  The extra normaliser element is, at the prime h,
// the translation by 1/h, i.e. the matrix (h 1; 0 h) of reduced norm h^2 and trace t; classically
// t = 2h and it is parabolic, but a Shimura curve with D > 1 is compact, so the global element is
// elliptic and its fixed points are CM points of discriminant t^2 - 4h^2 with |t| < 2h.  Those
// discriminants are NOT of the form -q or -4q for q | DN', so such a point is fixed by no
// Atkin-Lehner involution and is absent from the elliptic/AL enumeration of branchA4.m.
//
// For an image of order 2 in the Galois group the element must square into Q^* O^*, i.e. t = 0,
// giving d = -4h^2: d = -36 for h = 3 and d = -16 for h = 2.  This script asks, for each
// obstructed base, whether the base carries star points of those discriminants and whether adding
// them closes Riemann-Hurwitz.
AttachSpec("ShimuraQuotients.spec");
hall := func<M | [d : d in Divisors(M) | GCD(d, M div d) eq 1]>;
for b in [[10,9],[14,9],[22,9],[15,8],[15,4],[21,4],[33,4]] do
    D := b[1]; N := b[2];
    h := (N mod 9 eq 0) select 3 else 2;
    Np := N div h^2;
    W := hall(D*Np);
    psi := func<M | M * &*[Rationals() | 1 + 1/q : q in PrimeDivisors(M)]>;
    degX := Integers() ! (psi(N) / psi(Np) * #W / #hall(D*N));   // index of Gamma_0(N) in Gamma_0(N'), then the two star quotients
    gtop := GenusShimuraCurveQuotient(D, N, Set(hall(D*N)));
    gbot := GenusShimuraCurveQuotient(D, Np, Set(hall(D*Np)));
    R := 2*gtop - 2 - degX*(2*gbot - 2);
    // candidate extra discriminants: t^2 - 4h^2 for |t| < 2h, discarding non-discriminants
    cands := [t^2 - 4*h^2 : t in [0..2*h-1]];
    cands := [d : d in cands | d lt 0 and (d mod 4 in [0,1])];
    printf "\n== %o_%o: degree %o over %o_%o, needs R = %o ; extra-element discriminants %o\n",
           D, N, degX, D, Np, R, cands;
    for d in cands do
        Rd := QuadraticOrder(BinaryQuadraticForms(d));
        n := NumberOfOptimalEmbeddings(Rd, D, Np);
        if n eq 0 then printf "   d = %-5o : no CM points on the base\n", d; continue; end if;
        fixers := [q : q in W | q ne 1 and IsDefined(NumFixedPointsByCMOrder(D, Np, q), d)
                                       and NumFixedPointsByCMOrder(D, Np, q)[d] eq n];
        fixv := [n] cat [IsDefined(NumFixedPointsByCMOrder(D, Np, q), d)
                         select NumFixedPointsByCMOrder(D, Np, q)[d] else 0 : q in W | q ne 1];
        printf "   d = %-5o : %o point(s) on the base, %o star point(s), fixed by %o\n",
               d, n, &+fixv / #W, fixers;
    end for;
end for;
