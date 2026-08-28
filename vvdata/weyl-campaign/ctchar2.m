// LEMMA half 2, the convention-free version.
//
// ctchar.m compared a_0(theta|gamma) against bare Kronecker symbols and matched
// nothing -- naturally: in half-integral weight the measured quantity carries
// the theta factor eps_d^{-3}(c|d) too, and ctTheta works in the DUAL (QSIGN=-1)
// normalisation besides.  Recalling the classical formula through two layers of
// convention is exactly how this campaign has burned time before.
//
// So ask the question that needs no convention at all.  Write
// nu_L(gamma) := a_0(theta_L | gamma) for gamma in Gamma_0(4DN); since
// a_0(theta_L) = 1 this IS the multiplier.  For the alternating sum of
// Theorem 9.11 to be a modular form at all, every lattice in it must share one
// multiplier.  Test that directly -- same code path, same gamma, ratios only:
//
//   (i)  is nu_L the same for all (Delta,R) in a base's combination?
//   (ii) is it the same ACROSS bases (it should be, since det A = 2(Delta R)^2
//        makes 2 det A a perfect square, so no dependence on Delta or R can
//        enter through the discriminant)?
//
// Where the pool test already succeeded (15_2, 33_2, 35_2, 10_7, 10_1, 14_1,
// 26_1) the common value IS the panel's conjugate character; (i) + (ii) then
// carry that to the bases whose pools are too small to decide.
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/fastcosets.m";
load "vvdata/weyl-campaign/etared.m";
load "vvdata/weyl-campaign/ctlat.m";
PREC := 60; CC := ComplexField(PREC); QSIGN := -1;
R8 := ComplexField(8);

terms := function(D, N)
    out := {@ @};
    for s in PrimeDivisors(D) do
        Dp := D div s;
        for T in Subsets(SequenceToSet(PrimeDivisors(N))) do
            key := &*[ Integers() | p : p in T ];
            Rr := N div key;
            if IsOdd(#T) then Include(~out, <Dp*s*key, Rr>);
            else Include(~out, <Dp*key, Rr>); Include(~out, <Dp*key, Rr*s>);
            end if;
        end for;
    end for;
    return [ t : t in out | #PrimeDivisors(t[1]) mod 2 eq 1 ];
end function;

bases := [ <15,2>, <21,2>, <22,3>, <33,2>, <35,2>, <55,2>, <10,7>, <34,3>,
           <10,1>, <14,1>, <22,1>, <26,1> ];
ds := [1, 3, 7, 9, 11, 13, 17, 19, 23];

NU := AssociativeArray();          // (base, d) -> the common value, if any
allsame := true;
for b in bases do
    D := b[1]; N := b[2]; M := 4*D*N;
    tl := terms(D, N);
    printf "\n=== %o_%o  4DN = %o  lattices %o ===\n", D, N, M,
        [ Sprintf("Gr(%o,%o)", t[1], t[2]) : t in tl ];
    // gamma = [a b; M d] in Gamma_0(M) <= Gamma_0(level(L)) for every L
    gams := [ ]; dsu := [ ];
    for d in ds do
        if GCD(d, M) ne 1 then continue; end if;
        a := InverseMod(d mod M, M);
        b2 := (a*d - 1) div M;
        gam := Matrix(Integers(), 2, 2, [a, b2, M, d]);
        if Determinant(gam) ne 1 then continue; end if;
        Append(~gams, gam); Append(~dsu, d);
    end for;
    words := [ VVSTWord(g) : g in gams ];
    vals := [ ];
    for t in tl do
        A := ChangeRing(grossGram(t[1], t[2]), Rationals());
        v := ctTheta(A, words, CC, QSIGN);
        Append(~vals, v);
        printf "  Gr(%-3o,%-3o)  %o\n", t[1], t[2],
            [ R8!v[i] : i in [1..Minimum(6,#v)] ];
    end for;
    same := true;
    for i in [2..#vals] do
        for j in [1..#words] do
            if Abs(vals[i][j] - vals[1][j]) gt 10^(-30) then same := false; end if;
        end for;
    end for;
    if not same then allsame := false; end if;
    printf "  --> all lattices share one multiplier: %o\n", same;
    if same then
        for j in [1..#words] do NU[<D*100+N, dsu[j]>] := vals[1][j]; end for;
    end if;
end for;

printf "\n=== is the multiplier the same ACROSS bases? (by d) ===\n";
for d in ds do
    ks := [ k : k in Keys(NU) | k[2] eq d ];
    if #ks le 1 then continue; end if;
    v0 := NU[ks[1]];
    agree := forall{ k : k in ks | Abs(NU[k] - v0) le 10^(-30) };
    printf "  d = %-3o over %o bases: %o   value %o\n", d, #ks,
        agree select "SAME" else "DIFFERS", R8!v0;
end for;
printf "\nevery base internally consistent: %o\n", allsame;
printf "DONE\n";
quit;
