// The kappa closed form: the numerically-pinned constant of cusp1/cusp3 IS the
// Dedekind-eta multiplier system, made exact.
//
// For monomial m_r = prod_d eta(d tau)^{r_d} (weight 1/2: sum r_d = 1) and coset
// word w with matrix g, [d 0;0 1] g = g_d [a_d b_d; 0 e_d]:
//   eta(d g tau) = eps(g_d) * ((c tau + d')/e_d * (-i))^{1/2}-bookkeeping * eta(sigma_d tau)
// (Apostol Thm 3.4), whence
//   kappa_{r,w} = zeta_w * prod_d eps(g_d)^{r_d} * e_d^{-r_d/2}
// with zeta_w INDEPENDENT of r.  Test: for every coset and every monomial pair,
//   kappa_{r,w}/kappa_{r',w} = prod_d eps(g_d)^{r_d - r'_d} e_d^{-(r_d - r'_d)/2}.
// Also dumps zeta_w (= kappa/closed-form) per coset for the rho-side correlation.
//
//   magma -b DD:=15 NN:=2 KEYS:=-2,-1,9,10,11,12,13,14,15 cusp4.m
AttachSpec("ShimuraQuotients.spec");

D := 15; N := 2; keystr := "-2,-1,9,10,11,12,13,14,15";
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;
if assigned KEYS then keystr := KEYS; end if;
keys := [StringToInteger(s) : s in Split(keystr, ",")];

PREC := 80;
CC := ComplexField(PREC); ii := CC.1; pi := Pi(CC);
ee := func< z | Exp(2*pi*ii*z) >;

Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
R := Parent(fs[keys[1]]); ds := R`ds;
M := IsOdd(D*N) select 4*D*N else 2*D*N;
printf "BASE %o %o  M = %o  ds = %o\n", D, N, M, ds;

reps := VVCosetReps(M);
words := [ VVSTWord(g) : g in reps ];
printf "#cosets = %o\n", #reps;

triang := function(g, d)
    g2 := Matrix(Integers(), 2, 2, [d*g[1][1], d*g[1][2], g[2][1], g[2][2]]);
    c1 := g2[1][1]; c2 := g2[2][1];
    h := GCD(c1, c2);
    p1 := c1 div h; p2 := c2 div h;
    gg, u, v := XGCD(p1, p2);
    error if gg ne 1, "not primitive";
    gd := Matrix(Integers(), 2, 2, [p1, -v, p2, u]);
    sd := gd^(-1) * g2;
    a := sd[1][1]; b := sd[1][2]; e := sd[2][2];
    if a lt 0 then a := -a; b := -b; e := -e; end if;
    if e lt 0 then gd := -gd; sd := -sd; a := sd[1][1]; b := sd[1][2]; e := sd[2][2]; end if;
    return a, b, e, gd;
end function;

slashdata := function(word, tau)
    z := tau; factor := CC!1;
    for i := #word to 1 by -1 do
        if word[i][1] eq "S" then factor /:= Sqrt(z); z := -1/z;
        else z := z + word[i][2]; end if;
    end for;
    return factor, z;
end function;

// Dedekind sum s(h, k) = sum_{i=1}^{k-1} ((i/k)) ((h i/k)), exact rational
dedsum := function(h, k)
    s := Rationals()!0;
    for i := 1 to k - 1 do
        x := Rationals()!i/k;
        y := Rationals()!(h*i)/k;
        y := y - Floor(y);
        if y ne 0 then s +:= (x - 1/2)*(y - 1/2); end if;
    end for;
    return s;
end function;

// --- the eta multiplier eps(gamma), Apostol Thm 3.4 convention:
//     for c > 0: eta(g tau) = eps(g) * (-i (c tau + d))^(1/2) * eta(tau)
//     eps(g) = exp(pi i ((a + d)/(12 c) + s(-d, c)))
//     c = 0, a = d = 1: eta(tau + b) = e(b/24) eta(tau);  a = d = -1: same (g ~ -g).
epseta := function(g)
    a := g[1][1]; b := g[1][2]; c := g[2][1]; d := g[2][2];
    if c lt 0 or (c eq 0 and d lt 0) then a := -a; b := -b; c := -c; d := -d; end if;
    if c eq 0 then return ee(CC!b/24), 0;  end if; // (eps, cpos)
    s := dedsum(-d, c);
    return Exp(pi*ii*( CC!(a + d)/(12*c) + CC!s )), 1;
end function;

// numeric SELF-TEST of the transformation formula: returns err
selftest := function(gd)
    z := CC!0.123 + CC!0.917*ii;
    a := gd[1][1]; b := gd[1][2]; c := gd[2][1]; d := gd[2][2];
    gz := (a*z + b)/(c*z + d);
    eps, cpos := epseta(gd);
    if cpos eq 0 then
        pred := eps * DedekindEta(z);
    else
        cc := c; dd := d;
        if c lt 0 or (c eq 0 and d lt 0) then cc := -c; dd := -d; end if;
        pred := eps * Sqrt(-ii*(cc*z + dd)) * DedekindEta(z);
    end if;
    return Abs(pred - DedekindEta(gz));
end function;

monos := {@ @};
for k in keys do for r in Exponents(fs[k]) do Include(~monos, r); end for; end for;
printf "#monomials = %o\n", #monos;

tau0 := CC!0.31 + CC!1.31*ii;
tau1 := CC!(-0.57) + CC!1.73*ii;

tested := {}; stfails := 0;
npass := 0; nfail := 0; nnu := 0; nufails := 0;
betacases := AssociativeArray();
for wi->w in words do
    g := VVWordMatrix(w);
    fac0, z0 := slashdata(w, tau0);
    fac1, z1 := slashdata(w, tau1);
    x0 := CC!g[2][1]*tau0 + CC!g[2][2];   // c tau + d' of the ORIGINAL coset matrix
    x1 := CC!g[2][1]*tau1 + CC!g[2][2];
    tri := [ ]; us := [ CC | ];
    for d in ds do
        a, b, e, gd := triang(g, d);
        key := <gd[1][1], gd[1][2], gd[2][1], gd[2][2]>;
        if key notin tested then
            Include(~tested, key);
            err := selftest(gd);
            if err gt 10^(-60) then
                stfails +:= 1;
                printf "SELFTEST FAIL %o  err %o\n", key, RealField(6)!err;
            end if;
        end if;
        // nu_d = eta(d g tau) / [ (-i(c tau + d'))^{1/2} e^{-1/2} eta(sigma_d tau) ]
        sig0 := (CC!a*tau0 + CC!b)/CC!e;
        sig1 := (CC!a*tau1 + CC!b)/CC!e;
        nu0 := DedekindEta(CC!d*z0) / ( Sqrt(-ii*x0) * CC!e^(-CC!1/2) * DedekindEta(sig0) );
        nu1 := DedekindEta(CC!d*z1) / ( Sqrt(-ii*x1) * CC!e^(-CC!1/2) * DedekindEta(sig1) );
        nnu +:= 1;
        if Abs(nu0 - nu1) gt 10^(-55) then
            nufails +:= 1;
            printf "NU NOT CONSTANT coset %o d=%o dev %o\n", wi, d,
                RealField(6)!Abs(nu0 - nu1);
        end if;
        // beta case table: nu / eps(gd normalized) should be a 4th root of unity
        eps, _ := epseta(gd);
        beta := nu0/eps;
        b4 := Round(4*Arg(beta)/(2*Pi(RealField(30))));
        ckey := <Sign(gd[2][1]), Sign(gd[2][2]), Sign(gd[1][1]),
                 Abs(beta^4 - 1) lt 10^(-50) select b4 mod 4 else -99>;
        if IsDefined(betacases, ckey) then betacases[ckey] +:= 1;
        else betacases[ckey] := 1; end if;
        Append(~tri, <a, b, e, gd>);
        Append(~us, nu0 * ee(CC!b/(24*e)) * CC!e^(-CC!1/2));
    end for;
    W := LCM([ tri[i][3] : i in [1..#ds] ]);
    // numeric kappa per monomial (as in cusp1: num/sfun at tau0)
    kaps := [ CC | ]; preds := [ CC | ];
    for ri->r in monos do
        num0 := fac0 * &*[ CC | DedekindEta(d*z0)^(r[i]) : i->d in ds | r[i] ne 0 ];
        L := &+[ Integers() | r[i]*tri[i][1]*(W div tri[i][3]) : i in [1..#ds] ];
        sf := ee(tau0*L/(24*W)) *
            &*[ CC | ( DedekindEta((tri[i][1]*tau0 + tri[i][2])/tri[i][3]) *
                       ee(-(tri[i][1]*tau0 + tri[i][2])/(24*tri[i][3])) )^(r[i])
                : i in [1..#ds] | r[i] ne 0 ];
        Append(~kaps, num0/sf);
        Append(~preds, &*[ CC | us[i]^(r[i]) : i in [1..#ds] | r[i] ne 0 ]);
    end for;
    // kappa = zeta_w * prod_d u_d^{r_d} with zeta_w = fac0 * (-i x0)^{1/2}-free:
    zet := kaps[1]/preds[1];
    ok := true;
    for ri := 2 to #monos do
        if Abs(kaps[ri]/preds[ri] - zet) gt 10^(-50) then
            ok := false;
            printf "RATIO FAIL coset %o mono %o: dev %o\n", wi, ri,
                RealField(6)!Abs(kaps[ri]/preds[ri] - zet);
        end if;
    end for;
    if ok then npass +:= 1; else nfail +:= 1; end if;
    // zeta_w should equal fac0*Sqrt(-i x0) up to nothing: check directly
    zpred := fac0 * Sqrt(-ii*x0);
    cc := g[2][1] mod M; dd := g[2][2] mod M;
    printf "ZETA c=%o d=%o gcd=%o  zeta=(%o,%o)  zeta/pred=(%o,%o)  ok=%o\n",
        cc, dd, GCD(cc, M),
        RealField(10)!Re(zet), RealField(10)!Im(zet),
        RealField(10)!Re(zet/zpred), RealField(10)!Im(zet/zpred), ok;
end for;
printf "COSETS: kappa = zeta_w prod u_d^{r_d} holds on %o / %o (selftest fails %o; nu pairs %o, nonconstant %o)\n",
    npass, npass + nfail, stfails, nnu, nufails;
printf "BETA CASES <sign c, sign d, sign a, beta as i^k (or -99)> -> count:\n";
for ckey in Sort([k : k in Keys(betacases)]) do
    printf "  %o -> %o\n", ckey, betacases[ckey];
end for;
quit;
