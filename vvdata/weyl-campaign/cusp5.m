// Per-MONOMIAL class sums: the derivation lever the kappa closed form unlocks.
// c_{eta*}(0) is linear in f, so evaluate the cusp-class sums T_g on each eta
// monomial m_r separately (86 on 15_2, vs 9 panel forms): together they pin the
// FULL linear functional pp |-> c_{eta*}(0), from which the two channels
// (N|m support, square-class weights, psi-twist, w_sq) are read off by
// restriction -- no aliasing left.
//
// Output: MONO lines (index, exponent vector), then per monomial the class sums
//   TG <mono index> g=<class> re im   (50 digits)
//
//   magma -b DD:=15 NN:=2 KEYS:=-2,-1,9,10,11,12,13,14,15 cusp5.m
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

Ld := ShimuraCurveLattice(D, N);
fftdata := VVWeilFFT(Ld, CC : Dual := true);
elts := fftdata[7]; i0 := fftdata[8];
Qr := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
qofidx := function(i)
    v := ChangeRing(elts[i]@@Ld`to_disc, Rationals());
    r := (v*Qr, v)/(2*dn^2); return r - Floor(r);
end function;
isoidx := [ i : i in [1..#elts] | qofidx(i) eq 0 ];
printf "#isotropic = %o (2N-1 = %o)\n", #isoidx, 2*N-1;

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

tau0 := CC!0.31 + CC!1.31*ii;
tau1 := CC!(-0.57) + CC!1.73*ii;

monos := {@ @};
for k in keys do for r in Exponents(fs[k]) do Include(~monos, r); end for; end for;
printf "#monomials = %o\n", #monos;
for ri->r in monos do
    printf "MONO %o %o\n", ri, [r[i] : i in [1..#ds]];
end for;

SS := PowerSeriesRing(CC); t := SS.1;
a0tab := [ [ CC!0 : r in monos ] : w in words ];
for wi->w in words do
    g := VVWordMatrix(w);
    tri := [ ];
    for d in ds do
        a, b, e, gd := triang(g, d);
        Append(~tri, <a, b, e, gd>);
    end for;
    W := LCM([ tri[i][3] : i in [1..#ds] ]);
    leads := [ &+[ Integers() | r[i]*tri[i][1]*(W div tri[i][3]) : i in [1..#ds] ]
               : r in monos ];
    depth := Maximum([ 0 ] cat [ -L : L in leads ]) + 1;
    units := [ ];
    for i->d in ds do
        a, b, e := Explode([tri[i][1], tri[i][2], tri[i][3]]);
        step := 24*a*(W div e);
        u := SS!1 + O(t^depth);
        n := 1;
        while n*step lt depth do
            u *:= 1 - ee(CC!(n*b/e))*t^(n*step);
            n +:= 1;
        end while;
        Append(~units, u);
    end for;
    fac0, z0 := slashdata(w, tau0);
    fac1, z1 := slashdata(w, tau1);
    for ri->r in monos do
        L := leads[ri];
        if L gt 0 then continue; end if;
        produ := &*[ SS | units[i]^(r[i]) : i in [1..#ds] | r[i] ne 0 ];
        c0 := Coefficient(produ, -L);
        if c0 eq 0 and -L gt 0 then a0tab[wi][ri] := CC!0; continue; end if;
        num0 := fac0 * &*[ CC | DedekindEta(d*z0)^(r[i]) : i->d in ds | r[i] ne 0 ];
        num1 := fac1 * &*[ CC | DedekindEta(d*z1)^(r[i]) : i->d in ds | r[i] ne 0 ];
        sfun := func< tau | ee(tau*L/(24*W)) *
            &*[ CC | ( DedekindEta((tri[i][1]*tau + tri[i][2])/tri[i][3]) *
                       ee(-(tri[i][1]*tau + tri[i][2])/(24*tri[i][3])) )^(r[i])
                : i in [1..#ds] | r[i] ne 0 ] >;
        k0 := num0 / sfun(tau0);
        k1 := num1 / sfun(tau1);
        error if Abs(k0 - k1) gt 10^(-30), "kappa not constant", wi, ri, Abs(k0-k1);
        a0tab[wi][ri] := k0 * c0;
    end for;
end for;
printf "a0 table built\n";

// per-monomial class sums over eta* = isoidx[2]
est := isoidx[2];
Tg := AssociativeArray();
for wi->w in words do
    if forall{ri : ri in [1..#monos] | Abs(a0tab[wi][ri]) lt 10^(-PREC+10)} then
        continue;
    end if;
    g := VVWordMatrix(w);
    rv := VVRhoInvE0FFT(fftdata, w);
    cls := GCD(g[2][1] mod M, M);
    for ri := 1 to #monos do
        if Abs(a0tab[wi][ri]) lt 10^(-PREC+10) then continue; end if;
        key := <cls, ri>;
        if IsDefined(Tg, key) then Tg[key] +:= rv[est] * a0tab[wi][ri];
        else Tg[key] := rv[est] * a0tab[wi][ri]; end if;
    end for;
end for;
R50 := RealField(50);
for key in Sort([k : k in Keys(Tg)]) do
    printf "TG %o g=%o %o %o\n", key[2], key[1], R50!Re(Tg[key]), R50!Im(Tg[key]);
end for;
printf "DONE\n";
quit;
