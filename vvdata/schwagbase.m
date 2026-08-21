// THE CLOSED-FORM TEST: Schwagenscheidt's oldform relation reduces the seed-eta* Eisenstein series
// to the SEED-0 series of the enlarged lattice L1 = <L, eta*>, giving
//     mult_f = -(1/(4(N-1))) * Sum_{n>0} Sum_{gamma in eta*perp} c_gamma(-n) * b^{B,0}_{gamma+H}(n)
// with B = L1'/L1 and b^{B,0} from Bruinier-Kuss Theorem 7 (rank 3, weight 3/2, dual rep).
//
// Bruinier-Kuss, m = rank = 3 (odd), k = 3/2, so 1 - m/2 - k = -2 and the local variable is p^{-2};
// the coefficient of q^n e_gbar (n > 0, n = -Q(gbar) mod 1  [dual rep]) is
//     b(gbar, n) = C_glob * n^{k-1} * L(chi_D0, k-1/2)/L(chi^2, 2k-1) *
//                  prod_{p | 2 Ng^2 n det} [ (1 - chi_D0(p) p^{1/2-k}) / (1 - p^{1-2k}) * Lp(p^{-2}) ]
//     Lp(X) = N(p^w) X^w + (1 - p^{m-1} X) Sum_{nu<w} N(p^nu) X^nu,   w = 1 + 2 v_p(2 Ng^2 n),
//     N(a)  = #{ r in L1/aL1 : Q(r - gbar) + n = 0 mod a },
//     D = 2 (-1)^((m+1)/2) Ng^2 n det(L1),  D0 = fundamental discriminant of D.
// With k = 3/2: L(chi_D0, 1) = 2 pi h(D0) / (w(D0) sqrt|D0|)  and  L(chi^2, 2) = zeta(2) = pi^2/6,
// so  b = C * sqrt(n / |D0|) * (h/w) * prod_p [ (1 - chi_D0(p)/p) / (1 - 1/p^2) * Lp(p^-2) ]
// up to one global constant C -- which the RATIO test does not need.
//
// SIGN/REP CONVENTION: the input F lives in the dual rep with the code's Gram negQ; the same Gram is
// used for L1 throughout, exactly as the pipeline does.

D := 15; N := 2;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;
Ld := ShimuraCurveLattice(D, N);
negQ := -ChangeRing(Ld`Q, Integers());
Qr := ChangeRing(Ld`Q, Rationals());
negQr := ChangeRing(negQ, Rationals());
denom := Ld`denom; to_disc := Ld`to_disc;
A := Ld`disc_grp;
M := IsOdd(D*N) select 4*D*N else 2*D*N;
zero := Vector(Rationals(), [0,0,0]);

// oracle / consistency-relation ground truth per base
truthtab := AssociativeArray();
truthtab[[15,2]]  := [<-2,2>, <-1,4>, <9,0>, <10,0>, <11,4>, <12,2>, <13,4>, <14,-2>, <15,2>];
truthtab[[6,5]]   := [<-2,0>, <-1,0>, <9,6>, <10,3>, <11,12>, <12,9>, <13,3>, <14,0>, <15,3>];
truthtab[[10,3]]  := [<-2,0>, <-1,0>, <9,6>, <10,3>, <11,3>, <12,3>, <13,3>, <14,0>, <15,6>];
truthtab[[14,5]]  := [<-2,1>, <-1,1>];
truth := AssociativeArray();
for t in truthtab[[D,N]] do truth[t[1]] := t[2]; end for;

// ---- eta* and the enlarged lattice L1 = <Z^3, eta*> --------------------------------------------
iso := [ eta : eta in A | not IsZero(eta) and
         IsIntegral((v*Qr, v)/(2*denom^2)) where v := ChangeRing(eta@@to_disc, Rationals()) ];
wstar := ChangeRing(iso[1]@@to_disc, Rationals())/denom;
Nb := Lcm([Denominator(x) : x in Eltseq(wstar)] cat [1]);   // order of eta* in A
printf "eta* = %o, order N_beta = %o (level N = %o)\n", wstar, Nb, N;
error if Nb ne N, "expected the isotropic coset to have order N";

// L1 basis: rows of B1 in ambient Q^3 coordinates (which carry Gram negQ)
L1big := VerticalJoin(IdentityMatrix(Rationals(), 3), Matrix(Rationals(), 1, 3, Eltseq(wstar)));
H1 := HermiteForm(ChangeRing(L1big * Nb, Integers()));
B1 := ChangeRing(Matrix(Integers(), [Eltseq(H1[i]) : i in [1..3]]), Rationals()) / Nb;
G1 := B1 * negQr * Transpose(B1);
error if not &and[ IsIntegral(G1[i][j]) : i, j in [1..3] ] or not &and[ IsEven(Integers()!G1[i][i]) : i in [1..3] ],
      "L1 is not even integral";
G1z := ChangeRing(G1, Integers());
det1 := Determinant(G1z);
printf "L1 built: det = %o (det L = %o, ratio %o)\n", det1, Determinant(negQ), Determinant(negQ)/det1;

// ---- the discriminant group of L1, realized concretely via the dual ----------------------------
// classes gamma+H for gamma in eta*perp: we only ever need, per pairing coset gamma (ambient vector
// w), its class in B = L1'/L1 -- concretely just the ambient vector w regarded mod L1 -- and its
// order Ng in B, plus counting numbers N_{w,n}(p^nu) on L1.  All computable from B1 and G1.
tobasis := B1^-1;   // ambient -> L1-basis coordinates

// order of w mod L1: smallest t with t*w in L1
orderB := function(w)
    c := w * tobasis;
    return Lcm([Denominator(x) : x in Eltseq(c)] cat [1]);
end function;

// KEY SIMPLIFICATION: at rank m = 3, weight k = 3/2, the Bruinier-Kuss local polynomial is
// evaluated at X = p^{1-m/2-k} = p^{-2}, where the factor (1 - p^{m-1}X) VANISHES.  So
//     L^{(p)}_{gamma,n}(p^{-2}) = N_{gamma,n}(p^{w}) / p^{2w}  =  the stabilized LOCAL DENSITY
// on gamma + L1 -- exactly the object TwoCosetW (validated 126/126 against the pipeline) computes.
// No w_p bookkeeping survives; the density is computed at a stabilized level and CHECKED stable.
densityL1 := function(w, n, p)
    c := w * tobasis;                                    // coset in L1-basis coordinates
    zz := Vector(Rationals(), [0,0,0]);
    // adaptive level: accept when two consecutive levels agree (densities are eventually constant;
    // starting low keeps p = 3 at 3^9-3^12 loops in the common case and escalates only when needed)
    k := p eq 2 select 5 else (p eq 3 select 3 else 2);
    prev := TwoCosetW(n, p, k, c, zz, G1z);
    for it in [1..4] do
        cur := TwoCosetW(n, p, k+it, c, zz, G1z);
        if cur eq prev then return Rationals()!cur; end if;
        prev := cur;
    end for;
    error Sprintf("density not stabilized at p=%o n=%o", p, n);
end function;

// b^{B,0}_{w}(n) up to ONE global constant (the ratio test needs no more)
bB0 := function(w, n)
    Ng := orderB(w);
    Dd := Integers()!(2 * Ng^2 * Numerator(n)*Denominator(n) * det1);   // up to squares: 2 Ng^2 n det
    D0 := FundamentalDiscriminant(Dd);
    error if D0 ge 0, "expected an imaginary quadratic character here";
    chi := KroneckerCharacter(D0);
    hh := ClassNumber(QuadraticField(D0));
    ww := #TorsionSubgroup(UnitGroup(QuadraticField(D0)));
    ok, s := IsSquare(Rationals()!(n/AbsoluteValue(D0)));
    error if not ok, "sqrt(n/|D0|) irrational -- discriminant bookkeeping is off";
    val := s * (Rationals()!hh/ww);
    for p in PrimeDivisors(2 * Ng^2 * Numerator(n) * det1) do
        val *:= (1 - Evaluate(chi, p)/p) / (1 - 1/(Rationals()!p)^2) * densityL1(w, n, p);
    end for;
    return val;
end function;

// ---- the input forms' principal parts, restricted to eta*perp ----------------------------------
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
ks := Sort([k : k in Keys(fs)]);

mod_M_to_vecs := AssociativeArray([0..M-1]);
for j in [0..M-1] do mod_M_to_vecs[j] := []; end for;
for eta in A do
    v := ChangeRing(eta@@to_disc, Rationals());
    res := M*((v*Qr, v)/(2*denom^2));
    if not IsIntegral(res) then continue; end if;
    Append(~mod_M_to_vecs[Integers()!res mod M], eta);
end for;

// pairing gamma must satisfy (gamma, eta*) in Z  (i.e. gamma in eta*perp)
inperp := func<w | IsIntegral((w * negQr, wstar))>;

printf "\n%-6o %-9o %-18o\n", "form", "truth", "candidate sum S_f";
res := AssociativeArray();
for kf in ks do
    foo := qExpansionAtoo(fs[kf], 1);
    f0  := qExpansionAt0(fs[kf], 1);
    S := Rationals()!0;
    for n in [1..-Valuation(foo)] do
        c := Coefficient(foo, -n);
        if c eq 0 then continue; end if;
        key := <Eltseq(zero), Rationals()!n>;
        if not IsDefined(cacheB, key) then cacheB[key] := bB0(zero, Rationals()!n); end if;
        S +:= c * cacheB[key];
    end for;
    for j in [1..-Valuation(f0)] do
        c := Coefficient(f0, -j);
        if c eq 0 then continue; end if;
        r := (Rationals()!j)/M;
        for eta in mod_M_to_vecs[j mod M] do
            w := ChangeRing(eta@@to_disc, Rationals())/denom;
            if not inperp(w) then continue; end if;
            key := <Eltseq(w), r>;
            if not IsDefined(cacheB, key) then cacheB[key] := bB0(w, r); end if;
            S +:= c * cacheB[key];
        end for;
    end for;
    res[kf] := S;
    printf "%-6o %-9o %-18o\n", kf, truth[kf], S;
end for;
printf "\nratios vs form -1 (truth ratio: mult/4):\n";
for kf in ks do
    printf "  form %-4o truth %-7o candidate %o\n", kf, truth[kf]/truth[-1],
           res[-1] eq 0 select "n/a" else Sprint(res[kf]/res[-1]);
end for;

// ---- the -E_{A,0} piece, on the SAME footing: seed-0 densities on L itself ---------------------
densityA := function(w, n, p)
    zz := Vector(Rationals(), [0,0,0]);
    k := p eq 2 select 5 else (p eq 3 select 3 else 2);
    prev := TwoCosetW(n, p, k, w, zz, negQ);
    for it in [1..4] do
        cur := TwoCosetW(n, p, k+it, w, zz, negQ);
        if cur eq prev then return Rationals()!cur; end if;
        prev := cur;
    end for;
    error Sprintf("A-density not stabilized at p=%o n=%o", p, n);
end function;
detA := Determinant(negQ);
orderA := func<w | Lcm([Denominator(x) : x in Eltseq(w)] cat [1])>;
bA0 := function(w, n)
    Ng := orderA(w);
    Dd := Integers()!(2 * Ng^2 * Numerator(n)*Denominator(n) * detA);
    D0 := FundamentalDiscriminant(Dd);
    error if D0 ge 0, "expected imaginary";
    chi := KroneckerCharacter(D0);
    hh := ClassNumber(QuadraticField(D0));
    ww := #TorsionSubgroup(UnitGroup(QuadraticField(D0)));
    ok, s := IsSquare(Rationals()!(n/AbsoluteValue(D0)));
    error if not ok, "sqrt irrational (A)";
    val := s * (Rationals()!hh/ww);
    for p in PrimeDivisors(2 * Ng^2 * Numerator(n) * detA) do
        val *:= (1 - Evaluate(chi, p)/p) / (1 - 1/(Rationals()!p)^2) * densityA(w, n, p);
    end for;
    return val;
end function;

// mult-candidate = -K [ S_B - S_A/N ]; with the full Schwagenscheidt/BK constants at k = 3/2,
// K = (12 * 2^(7/2) / sqrt|A|) / (4(N-1)) up to the i^kappa phase -- for 15_2, |A| = 1800:
// 12*2^(7/2)/(30 sqrt 2) = 32/5, so mult = -(8/5)(S_B - S_A/2) up to sign.  Print raw pieces too.
cacheA := AssociativeArray();
// the pairing constant: K = 12*2^(7/2)/(sqrt|B| * 4(N-1)), with the S_A term carrying the
// relative factor sqrt|B|/sqrt|A| = 1/N (hence S_A/N below).  |B| = |A|/N^2 = 2D^2, so
// sqrt|B| = D sqrt 2 and K = 12*8/(D * 4(N-1)) = 24/(D(N-1)).  Check: 15_2 -> 24/15 = 8/5.
K := Rationals()!24 / (D*(N-1));
printf "
pairing constant K = 24/(DN(N-1)) = %o;  candidates are -K (S_B - S_A/N) and +...
", K;
printf "
%-6o %-8o %-12o %-12o %-16o %-16o
",
       "form", "truth", "S_B", "S_A", "-K(SB-SA/N)", "+K(SB-SA/N)";
for kf in ks do
    foo := qExpansionAtoo(fs[kf], 1);
    f0  := qExpansionAt0(fs[kf], 1);
    SA := Rationals()!0;
    for n in [1..-Valuation(foo)] do
        c := Coefficient(foo, -n);
        if c eq 0 then continue; end if;
        key := <Eltseq(zero), Rationals()!n>;
        if not IsDefined(cacheA, key) then cacheA[key] := bA0(zero, Rationals()!n); end if;
        SA +:= c * cacheA[key];
    end for;
    for j in [1..-Valuation(f0)] do
        c := Coefficient(f0, -j);
        if c eq 0 then continue; end if;
        r := (Rationals()!j)/M;
        for eta in mod_M_to_vecs[j mod M] do
            w := ChangeRing(eta@@to_disc, Rationals())/denom;
            key := <Eltseq(w), r>;
            if not IsDefined(cacheA, key) then cacheA[key] := bA0(w, r); end if;
            SA +:= c * cacheA[key];
        end for;
    end for;
    v := K*(res[kf] - SA/N);
    printf "%-6o %-8o %-12o %-12o %-16o %-16o
", kf, truth[kf], res[kf], SA, -v, v;
end for;
quit;
