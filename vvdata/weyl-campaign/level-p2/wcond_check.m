// Check lem:Wcond (the closed form of the m = 0 counting series at an odd conductor prime) against
// brute-force counting on the ACTUAL plane L_- = L cap lambda^perp of X_0^10(3), p = 3:
// d = -72, -648 (3 split in Q(sqrt -8), k = 1, 2) and -387, -603 (3 inert in Q(sqrt -43), Q(sqrt -67), k = 1).
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
D := 10; N := 3; p := 3; JMAX := 6;
Ld := ShimuraCurveLattice(D, N);
Qg := ChangeRing(Ld`Q, Rationals()); n := Nrows(Qg);
e := func< i | Vector(Rationals(), [j eq i select 1 else 0 : j in [1..n]]) >;
inLdual := func< v | forall{i : i in [1..n] | IsIntegral((v*Qg, e(i)))} >;
closed0 := func< j, k, eps | j le 2*k select p^(j div 2) else (1+eps)*(p-1)*p^(k-1)*Ceiling((j-2*k)/2) + p^(k - (j mod 2)) >;
closedr := func< j, rho, eps | j le 2*rho select p^(j div 2) else (1+eps)*p^rho >;
for d in [-72, -648, -387, -603] do
    lam := ChangeRing(ElementOfNorm(Ld`Q, -d, Ld`O, Ld`basis_L), Rationals());
    Mx := Matrix(Integers(), n, 1, [Integers() | (e(i)*Qg, lam) : i in [1..n]]);
    K := ChangeRing(KernelMatrix(Mx), Rationals());
    gram := K*Qg*Transpose(K);
    d0 := FundamentalDiscriminant(d); f := Isqrt(d div d0); k := Valuation(f, p); eps := KroneckerSymbol(d0, p);
    printf "\n=== d = %o  d0 %o  f %o  k %o  eps %o\n", d, d0, f, k, eps;
    A11 := Integers()!(gram[1,1]/2); A12 := Integers()!gram[1,2]; A22 := Integers()!(gram[2,2]/2);
    Sm, Pm := SmithForm(ChangeRing(gram, Integers())); Pq := ChangeRing(Pm, Rationals());
    s1 := Sm[1,1]; s2 := Sm[2,2];
    // p-primary integral cosets
    keys := {};
    for a in [0..s1-1], b in [0..s2-1] do
        v := Vector(Rationals(), [a/s1, b/s2]) * Pq;
        key := [c - Floor(c) : c in Eltseq(v)];
        if exists{c : c in key | Denominator(c) ne p^Valuation(Denominator(c), p)} then continue; end if;
        if not IsIntegral((Vector(Rationals(), key)*gram, Vector(Rationals(), key))/2) then continue; end if;
        Include(~keys, key);
    end for;
    ok := true;
    for key in keys do
        a := Maximum([Valuation(Denominator(c), p) : c in key] cat [0]);
        r1 := Integers()!(key[1]*p^a); r2 := Integers()!(key[2]*p^a);
        alphas := [];
        for j in [1..JMAX] do
            cnt := 0; md := p^(j + 2*a); sh := p^a;
            for u in [0..p^j-1] do U := r1 + sh*u;
                for t in [0..p^j-1] do T := r2 + sh*t;
                    if (A11*U^2 + A12*U*T + A22*T^2) mod md eq 0 then cnt +:= 1; end if;
                end for;
            end for;
            Append(~alphas, cnt / p^j);
        end for;
        nu := Vector(Rationals(), key) * K;
        if key eq [0,0] then
            pred := [closed0(j, k, eps) : j in [1..JMAX]];
            printf "  zero coset        counted %o  lemma %o  %o\n", alphas, pred, alphas eq pred select "OK" else "MISMATCH";
        else
            // depth rho: the coset has order p^(k - rho) in L_-^v/L_-
            ord := Minimum([m : m in [1..2*k] | forall{c : c in key | IsIntegral(c*p^m)}]);
            rho := k - ord;
            pred := [closedr(j, rho, eps) : j in [1..JMAX]];
            printf "  coset %-12o rho %o in L^v %-5o counted %o  lemma %o  %o\n", key, rho, inLdual(nu), alphas, pred, alphas eq pred select "OK" else "MISMATCH";
        end if;
        if alphas ne pred then ok := false; end if;
    end for;
    printf "  %o integral p-primary cosets (lemma: p^k = %o), all %o\n", #keys, p^k, ok select "OK" else "MISMATCH";
end for;
