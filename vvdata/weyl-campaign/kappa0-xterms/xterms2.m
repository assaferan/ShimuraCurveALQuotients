// Per (form, firing discriminant d): the order of psi_F at tau_d and the x != 0 part of the dropped
// m = 0 sum, both from the two cusp expansions via [GY, Lemma 24] and [Yang, Thm A(2)].
//   L_+ = Z lambda0, Q(lambda0) = |d| (d = 1 mod 4) or |d|/4 (d = 0 mod 4);
//   L^v cap Q lambda = Z lambda0/g, g = 2^e * prod_{p odd, p | (D,d)} p, e = 1 iff d = 1 mod 4 or 2 | D
//   (d = 0 mod 4 with 2 \nmid DN: e = 0); x = j lambda0/g, Q(x) = j^2 Q(lambda0)/g^2.
//   coefficient at [x]: cusp-0 coefficient at exponent -Q(x) (Lemma 24), plus the oo-coefficient
//   when x in L (g | j).  ord_{tau_d} psi = sum_{j >= 1} c_[x](-Q(x))  (Thm A(2), up to a common
//   positive factor); extra = sum_{x != 0} c_[x+nu](-Q(x)) = 2 * (cusp-0 parts only).
// A pair with ord = 0 (Theorem B applies) and extra != 0 is a DISCREPANCY with mult = c_eta(0)/2.
AttachSpec("ShimuraQuotients.spec");
D := StringToInteger(D_s); N := StringToInteger(N_s);
M := IsOdd(D*N) select 4*D*N else 2*D*N;
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
ks := Sort([k : k in Keys(fs)]);
printf "X_0^%o(%o), M = %o\n", D, N, M;
maxoo := Maximum([-Valuation(qExpansionAtoo(fs[k], 80)) : k in ks]);
// firing fundamental discriminants with CM(d) nonempty: N splits, p | D not split
ds := [d : d in [-4*maxoo-4..-3] | IsFundamentalDiscriminant(d) and d mod N ne 0
        and KroneckerSymbol(d, N) eq 1
        and forall{p : p in PrimeDivisors(D) | KroneckerSymbol(d, p) ne 1}];
ndisc := 0;
for k in ks do
    f := fs[k];
    foo := qExpansionAtoo(f, 80); f0 := qExpansionAt0(f, 80);
    v0 := -Valuation(f0); voo := -Valuation(foo);
    out := [];
    for d in ds do
        Q0 := (d mod 4 eq 0) select -d div 4 else -d;
        e := (d mod 4 eq 1 or IsEven(D)) select 1 else 0;
        g := 2^e * &*[Integers() | p : p in PrimeDivisors(D) | IsOdd(p) and d mod p eq 0];
        ord := 0; extra := 0;
        for j in [1..g*Floor(Sqrt(voo/Q0)) + g] do
            Qx := j^2 * Q0 / g^2;
            if Qx gt Maximum(voo, v0/M) then break; end if;
            ex := M*Qx;
            error if not IsIntegral(ex), "non-integral exponent", d, j;
            c0part := Coefficient(f0, -Integers()!ex);
            ord +:= c0part; extra +:= 2*c0part;
            if j mod g eq 0 and IsIntegral(Qx) then ord +:= Coefficient(foo, -Integers()!Qx); end if;
        end for;
        if extra ne 0 or ord ne 0 then
            Append(~out, <d, ord, extra>);
            if ord eq 0 and extra ne 0 then ndisc +:= 1; end if;
        end if;
    end for;
    printf "form %o (ord_oo %o, ord_0 %o/%o): (d, ord at tau_d, extra) = %o\n", k, voo, v0, M,
           [t : t in out | t[3] ne 0];
    printf "    on divisor only: %o\n", [t[1] : t in out | t[3] eq 0];
end for;
printf "DISCREPANCIES (ord = 0, extra != 0): %o\nXTERMS2_DONE\n", ndisc;
