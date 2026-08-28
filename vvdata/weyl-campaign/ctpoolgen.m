// LEMMA, half 2: the CHARACTER.  Does the Gross combination lie in the base's
// pool of holomorphic weight-3/2 eta quotients WITH THE PANEL CHARACTER?
//
// That pool is, per enum32.py, exactly the holomorphic same-character
// weight-3/2 eta quotients at the base's level.  So membership is a direct
// verification that the Gross combination carries the panel character -- i.e.
// that f * (Gross combination) has trivial multiplier, the hypothesis of the
// residue-pairing theorem.  Verified at 15_2 in ctpool15.log; this generalises
// it across every base with a banked pool.
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/ctlat.m";
DEPTH := 80;
R<q> := PowerSeriesRing(Rationals(), DEPTH);

massmul := function(Dp, Rl)
    ps := PrimeDivisors(Dp);
    m := Rationals()!(&*[ Integers() | p-1 : p in ps ]) / (48*2^(#ps-1));
    for p in PrimeDivisors(Rl) do m *:= Rationals()!(p+1)/2; end for;
    return m;
end function;

thetabar := function(Dp, Rl)
    A := ChangeRing(grossGram(Dp, Rl), Rationals());
    L := LatticeWithGram(A);
    reps := Representatives(Genus(L));
    raw := [ Rationals() | 0 : m in [0..2*DEPTH] ];
    den := Rationals()!0;
    for Lr in reps do
        w := Rationals()!1/#AutomorphismGroup(Lr);
        T := ThetaSeries(Lr, 2*DEPTH);
        for m in [0..2*DEPTH] do raw[m+1] +:= w*Coefficient(T, m); end for;
        den +:= w;
    end for;
    raw := [ x/den : x in raw ];
    error if exists{ m : m in [1..2*DEPTH by 2] | raw[m+1] ne 0 },
        "odd index nonzero -- norm convention wrong", Dp, Rl;
    return [ raw[2*m+1] : m in [0..DEPTH-1] ];
end function;

//         D,   N,  pool file suffix
bases := [ <15,2>, <21,2>, <22,3>, <33,2>, <35,2>, <55,2>, <77,2>, <10,7>,
           <34,3>, <10,1>, <14,1>, <22,1>, <26,1> ];

nok := 0; nrun := 0;
for b in bases do
    D := b[1]; N := b[2]; DN := D*N;
    M := IsOdd(DN) select 4*DN else 2*DN;
    ds := Divisors(M);
    fn := Sprintf("vvdata/weyl-campaign/epool_%o_%o.txt", D, N);
    ok, raw := ReadTest(fn);
    if not ok then printf "\n=== %o_%o: no pool file, skipped\n", D, N; continue; end if;
    Epool := [ [ Integers() | x : x in r ] : r in eval raw ];

    eta_unit := function(d0, e)
        ser := R!1; n := 1;
        while d0*n lt DEPTH do ser *:= (1 - q^(d0*n))^e; n +:= 1; end while;
        return ser;
    end function;
    mono := function(rr)
        s1 := &+[ ds[i]*rr[i] : i in [1..#ds] ];
        error if s1 mod 24 ne 0, "monomial not integral in q", rr;
        ser := R!1;
        for i in [1..#ds] do
            if rr[i] ne 0 then ser *:= eta_unit(ds[i], rr[i]); end if;
        end for;
        return q^(s1 div 24) * ser;
    end function;

    printf "\n=== %o_%o  M = %o  #pool = %o ===\n", D, N, M, #Epool;
    P := [ mono(r) : r in Epool ];
    A := Matrix(Rationals(), DEPTH, #P,
                [ [ Coefficient(P[j], m) : j in [1..#P] ] : m in [0..DEPTH-1] ]);

    for s in PrimeDivisors(D) do
        Dp := D div s;
        if #PrimeDivisors(Dp) mod 2 eq 0 then continue; end if;   // need omega odd
        // S = {} at N = 1, S = {N} at N prime: the two cases the campaign has
        m1 := massmul(Dp, N); m2 := massmul(Dp, N*s);
        t1 := thetabar(Dp, N); t2 := thetabar(Dp, N*s);
        Tel := [ Rationals() | (m2*t2[i] - m1*t1[i])/(m2 - m1) : i in [1..DEPTH] ];
        if N eq 1 then
            Gv := Tel;  desc := Sprintf("Tel_%o(%o;1)", s, Dp);
        else
            t3 := thetabar(DN, 1);
            Gv := [ Rationals() | (1/2)*Tel[i] - (1/2)*t3[i] : i in [1..DEPTH] ];
            desc := Sprintf("1/2 Tel_%o(%o;%o) - 1/2 tb(%o,1)", s, Dp, N, DN);
        end if;
        okc, sol := IsConsistent(Transpose(A), Vector(Rationals(), Gv));
        nrun +:= 1; if okc then nok +:= 1; end if;
        printf "  s=%-3o %-34o in pool span: %o", s, desc, okc;
        if okc then
            chk := [ &+[ sol[j]*Coefficient(P[j], m) : j in [1..#P] ] : m in [0..DEPTH-1] ];
            printf "   (exact: %o, %o nonzero coords)", chk eq Gv,
                #[ j : j in [1..#P] | sol[j] ne 0 ];
        end if;
        printf "\n";
    end for;
end for;
printf "\n%o of %o (base, s) instances lie in the panel-character pool\n", nok, nrun;
printf "DONE\n";
quit;
