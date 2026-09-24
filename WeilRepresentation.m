// The Weil representation rho_L of the discriminant form of an even lattice L, and (later) the
// vector-valued Borcherds input form F_f = sum_gamma f|_gamma rho_L(gamma^{-1}) e_0.
//
// Here L is the rank-3, signature (1,2) trace-zero lattice of an Eichler order of level N in the
// quaternion algebra of discriminant D (see QuaternionLatticeData.m). rho_L acts on the group
// algebra C[L^vee/L]; |L^vee/L| = det(Q) = 2*(D*N)^2. This file provides rho_L(S), rho_L(T) as
// matrices over a cyclotomic field, used as a verification oracle for the c_eta coefficients that
// enter Schofer's formula (in particular the outer m=0 constant term at level primes p|N).
//
// Conventions (Borcherds / Bruinier), signature (b+,b-) = (1,2):
//   rho_L(T) e_eta = e(nm(eta)) e_eta,           nm(eta) = <eta,eta>/2 = Q(eta)  (mod 1)
//   rho_L(S) e_eta = ( e((b- - b+)/8) / sqrt|L^vee/L| ) * sum_{eta'} e(-<eta,eta'>) e_{eta'}
// with e(x) = exp(2 pi i x). For (1,2): e((b- - b+)/8) = e(1/8), and sqrt|L^vee/L| = D*N*sqrt(2).

// The positive-real square root of a squarefree m > 0 inside K = Q(zeta_n0), via quadratic Gauss
// sums. (IsSquare returns an arbitrary Galois-conjugate root; the metaplectic Weil representation
// requires the positive real root -- the braid relation (ST)^3 = S^2 has an odd number of S factors
// and so is sensitive to the sign of 1/sqrt|L^vee/L|.) Requires 8 | n0 if 2|m, p | n0 for odd p|m,
// and 4 | n0 if some p == 3 mod 4 divides m.
function positive_real_sqrt(K, z, n0, m)
    ee := func<a | z^(Integers()!((a - Floor(a))*n0))>;   // e(a) = exp(2 pi i a)
    r := K!1; mm := m;
    if IsEven(mm) then
        r *:= ee(1/8) + ee(-1/8);          // zeta_8 + zeta_8^{-1} = sqrt(2) > 0
        mm := mm div 2;
    end if;
    for p in PrimeDivisors(mm) do          // mm odd squarefree
        g := &+[ KroneckerSymbol(a, p) * ee(a/p) : a in [1..p-1] ];   // Gauss sum
        // g = sqrt(p) if p == 1 (4), g = i*sqrt(p) if p == 3 (4)
        r *:= (p mod 4 eq 1) select g else g * ee(-1/4);
    end for;
    return r;
end function;

intrinsic WeilRepresentationST(Ld::QuaternionLatticeData) -> AlgMatElt, AlgMatElt, SeqEnum, FldCyc
{The Weil representation matrices rho_L(S) and rho_L(T) of the discriminant form (L^vee/L, nm) of
 the signature (1,2) lattice L of Ld, over a cyclotomic field. Also returns the ordered list of
 discriminant-group elements (indexing the rows/columns of the matrices, in the same order) and the
 cyclotomic field.}
    D := Ld`D; N := Ld`N;
    dg := Ld`disc_grp; Q := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
    // Level of the lattice: the smallest M with M*nm(eta) in Z and M*<eta,eta'> in Z for all eta.
    M := IsOdd(D*N) select 4*D*N else 2*D*N;
    // |L^vee/L| = det(Q); its 1/sqrt lives in Q(zeta_{4*sqfree}). Enlarge n0 to hold it (and the
    // 8 for the metaplectic S-phase e(1/8)). SquarefreeFactorization(n) = (squarefree part, sqrt of
    // square part), i.e. n = sqfree * sq^2.
    sqfree, sq := SquarefreeFactorization(Integers()!Determinant(Q));
    n0 := Lcm([M, 8, 4*sqfree]);
    K<z> := CyclotomicField(n0);
    ee := func<a | z^(Integers()!((a - Floor(a))*n0))>;    // e(a) = exp(2 pi i a), a rational

    elts := [g : g in dg]; n := #elts;
    vs := [ChangeRing(g@@Ld`to_disc, Rationals()) : g in elts];   // lifts to the scaled dual
    nm := [ (vs[i]*Q, vs[i])/(2*dn^2) : i in [1..n] ];            // nm(eta_i) = Q(eta_i) mod 1

    // 1/sqrt|L^vee/L| in K (|L^vee/L| = n = det(Q) = #dg = sqfree*sq^2), using the POSITIVE real root.
    invsqrt := 1/(sq * positive_real_sqrt(K, z, n0, sqfree));
    phase := ee(1/8);                              // e((b- - b+)/8), (b+,b-)=(1,2)

    T := DiagonalMatrix(K, [ee(nm[i]) : i in [1..n]]);
    S := ZeroMatrix(K, n, n);
    for i in [1..n] do
        for j in [1..n] do
            S[i][j] := invsqrt * phase * ee(-(vs[i]*Q, vs[j])/dn^2);
        end for;
    end for;
    return S, T, elts, K;
end intrinsic;
