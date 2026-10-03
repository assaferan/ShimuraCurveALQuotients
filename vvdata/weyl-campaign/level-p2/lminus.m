// Sachi's review: at d = -12 on X_0^15(2) the negative plane L_- is said to have only ONE nonzero
// isotropic coset, not the 2N-2 = 2 the proof gives when N does not divide d.  d = -12 has
// conductor 2, so N = 2 divides the CONDUCTOR while missing the fundamental discriminant -3 --
// exactly the case the code admits and the proof excludes.  Computed from the lattice itself:
// L_- = L cap lambda^perp, then the isotropic cosets of L_-^v / L_-.
AttachSpec("ShimuraQuotients.spec");
D := 15; N := 2;
Ld := ShimuraCurveLattice(D, N);
Q := ChangeRing(Ld`Q, Rationals());
e := func< i | Vector(Rationals(), [j eq i select 1 else 0 : j in [1..3]]) >;
for d in [-7, -15, -60, -12] do
    ok := true; lam := 0;
    try lam := ElementOfNorm(Ld`Q, -d, Ld`O, Ld`basis_L); catch err ok := false; end try;
    if not ok then printf "d = %-5o : no optimal embedding on this base\n", d; continue; end if;
    lv := ChangeRing(lam, Rationals());
    Mx := Matrix(Integers(), 3, 1, [Integers() | (e(i)*Q, lv) : i in [1..3]]);
    K := KernelMatrix(Mx);
    n := Nrows(K);
    gram := Matrix(Rationals(), n, n, [(ChangeRing(K[i], Rationals())*Q, ChangeRing(K[j], Rationals())) : i, j in [1..n]]);
    det := Integers() ! Determinant(gram);
    gi := gram^(-1);
    // dual vectors v with v*gram integral, modulo Z^n: enumerate v = w/|det|
    dd := Abs(det);
    seen := {}; iso := 0;
    for w in CartesianPower([0..dd-1], n) do
        v := Vector(Rationals(), [w[i]/dd : i in [1..n]]);
        if not forall{c : c in Eltseq(v*gram) | IsIntegral(c)} then continue; end if;
        key := [c - Floor(c) : c in Eltseq(v)];
        if key in seen then continue; end if;
        Include(~seen, key);
        r := (v*gram, v)/2;
        if r eq Floor(r) then iso +:= 1; end if;
    end for;
    fund := IsFundamentalDiscriminant(d);
    cond := Integers() ! Sqrt(d / FundamentalDiscriminant(d));
    printf "d = %-5o fundamental %-5o conductor %-3o 2-part of det %o : nonzero isotropic cosets of L_- = %o\n",
           d, fund, cond, Valuation(det, 2), iso - 1;
end for;
