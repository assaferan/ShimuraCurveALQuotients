// cusp6 = cusp5 + the two dumps that let the analysis run without aliasing:
//   FORMC: the panel forms' monomial coefficients (exact rationals);
//   PP: per cusp class (canonical coset) and monomial, the FULL principal part
//       kappa * Coeff(produ, j) at every pole exponent (L + j < 0), plus the
//       constant term (L + j = 0), 50 digits. Exponents reported as (L+j) in
//       the 1/(24 W_g) grid together with W_g, so analysis can match slots
//       across monomials.
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
//   magma -b DD:=15 NN:=2 KEYS:=-2,-1,9,10,11,12,13,14,15 cusp6.m
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
nPanel := #monos;
extraseqs := [ [-8,3,1,1,5,2,-3,0,-2,3,5,-6],[-8,3,2,0,3,2,1,0,-1,2,1,-4],[-8,3,2,0,6,2,1,0,-2,-1,1,-3],[-8,3,2,1,0,1,2,0,0,1,-1,0],[-8,3,2,1,6,1,0,0,-2,-1,1,-2],[-8,3,2,2,3,0,1,0,-1,-2,-1,2],[-8,4,2,-1,6,-2,-1,6,-2,1,4,-8],[-8,4,2,0,6,-3,-2,6,-2,1,4,-7],[-8,4,3,-1,4,0,-1,1,-1,3,2,-5],[-8,4,3,0,7,-1,-2,1,-2,0,2,-3],[-8,4,3,1,4,-2,-1,1,-1,-1,0,1],[-8,5,3,-1,7,-5,-4,7,-2,2,5,-8],[-8,5,4,-1,8,-3,-4,2,-2,1,3,-4],[-8,8,2,-2,3,-1,-1,0,-1,-1,-1,3],[-7,2,1,0,2,4,3,-1,0,1,0,-4],[-7,3,2,0,3,1,0,0,0,2,1,-4],[-7,3,2,1,0,0,1,0,1,1,-1,0],[-7,4,3,-1,4,-1,-2,1,0,3,2,-5],[-6,3,0,-1,3,1,0,2,3,2,1,-7],[-6,6,1,-2,-2,-2,3,2,5,0,-2,-2],[-6,7,1,-4,-2,-5,2,8,5,2,1,-8],[-5,1,0,0,2,3,5,0,-1,0,2,-6],[-5,3,0,-1,0,0,1,2,5,1,-1,-4],[-5,3,0,-1,3,0,-1,2,4,2,1,-7],[-5,3,1,-1,6,0,-3,1,-2,3,3,-5],[-5,4,1,-1,-3,0,8,0,1,-2,-1,-1],[-5,5,1,-3,-3,-3,7,6,1,0,2,-7],[-4,1,-1,-1,-2,3,6,3,1,-1,-1,-3],[-4,1,-1,-1,1,3,4,3,0,0,1,-6],[-4,1,0,0,-1,2,6,0,1,-1,0,-3],[-4,1,0,0,2,2,4,0,0,0,2,-6],[-4,1,0,1,-1,3,4,-2,1,-2,-3,3],[-4,2,0,-1,-1,0,3,4,1,0,0,-3],[-4,2,0,-1,2,2,-1,0,0,5,2,-6],[-4,2,1,0,0,1,2,-1,1,-1,-2,2],[-4,3,0,-1,3,-1,-2,2,5,2,1,-7],[-4,3,1,-2,0,-2,1,5,1,1,1,-4],[-4,3,1,-1,0,-1,-2,1,1,5,1,-3],[-4,3,1,-1,3,-1,-2,1,0,2,1,-2],[-4,3,1,-1,6,-1,-4,1,-1,3,3,-5],[-4,4,1,-3,3,-4,-3,7,0,4,4,-8],[-4,5,2,-2,5,-2,-3,2,-3,2,7,-8],[-3,0,-3,-1,2,8,-1,2,0,2,2,-7],[-3,1,-2,-1,0,5,-1,1,1,3,2,-5],[-3,1,0,1,3,1,-2,1,-2,1,5,-5],[-3,1,0,1,6,1,-4,1,-3,2,7,-8],[-3,1,1,0,1,-1,6,1,-1,-2,3,-5],[-3,2,1,-2,-2,-4,7,7,0,-1,4,-8],[-3,2,2,-1,5,-1,-2,2,-2,2,4,-7],[-3,4,2,-1,2,-2,1,1,-1,0,2,-4],[-3,5,2,-3,-1,-5,2,7,0,1,3,-7],[-2,0,0,0,0,3,4,0,0,-1,0,-3],[-2,0,0,0,3,3,2,0,-1,0,2,-6],[-2,0,0,1,0,2,3,0,0,-1,0,-2],[-2,0,0,1,3,2,1,0,-1,0,2,-5],[-2,0,0,2,3,1,0,0,-1,0,2,-4],[-2,1,0,-1,0,-1,2,6,0,1,3,-8],[-2,1,0,0,0,-2,1,6,0,1,3,-7],[-2,1,1,0,1,0,1,1,0,0,1,-3],[-2,1,1,0,4,0,-1,1,-1,1,3,-6],[-2,1,1,1,4,-1,-2,1,-1,1,3,-5],[-2,2,1,-1,1,-4,-1,7,0,2,4,-8],[-2,2,2,-1,2,-2,-1,2,0,1,2,-4],[-2,2,2,-1,5,-2,-3,2,-1,2,4,-7],[-2,6,0,-4,0,-3,-1,6,0,2,3,-6],[0,0,-2,-1,-3,2,3,2,5,2,0,-7],[1,-2,-2,0,-4,4,8,0,1,0,1,-6],[1,0,-1,-1,0,1,0,1,0,3,2,-5],[2,-1,-1,0,-3,0,4,1,2,1,2,-6],[2,0,-2,-1,-3,0,1,2,7,2,0,-7],[2,2,0,-2,-1,-1,0,2,-1,2,6,-8],[3,-1,0,-1,-1,0,1,2,0,2,3,-7] ];
U := Universe(monos);
for e in extraseqs do Include(~monos, U!e); end for;
printf "#monomials = %o (panel %o + extra %o)\n", #monos, nPanel, #monos - nPanel;
for ri->r in monos do
    printf "MONO %o %o\n", ri, [r[i] : i in [1..#ds]];
end for;
for k in keys do
    f := fs[k];
    printf "FORMC %o %o\n", k,
        [ <Index(monos, r), f`coeffs[r]> : r in Exponents(f) ];
end for;

SS := PowerSeriesRing(CC); t := SS.1;
a0tab := [ [ CC!0 : r in monos ] : w in words ];
ppdumped := {};
for wi->w in words do
    g := VVWordMatrix(w);
    cls0 := GCD(g[2][1] mod M, M);
    dopp := cls0 notin ppdumped;
    if dopp then Include(~ppdumped, cls0); end if;
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
        if c0 eq 0 and -L gt 0 and forall{j : j in [0..-L] | Coefficient(produ, j) eq 0} then
            a0tab[wi][ri] := CC!0; continue;
        end if;
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
        if dopp then
            R50 := RealField(50);
            for j := 0 to -L do
                cj := Coefficient(produ, j);
                if cj eq 0 then continue; end if;
                v := k0 * cj;
                printf "PP %o g=%o W=%o e=%o %o %o\n", ri, cls0, W, L + j,
                    R50!Re(v), R50!Im(v);
            end for;
        end if;
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
// per-form assembly self-check (same-run): c_eta*(0) per panel form
for k in keys do
    f := fs[k];
    tot := CC!0;
    for key in Keys(Tg) do
        if key[2] gt nPanel then continue; end if;
        r := monos[key[2]];
        if IsDefined(f`coeffs, r) then tot +:= f`coeffs[r] * Tg[key]; end if;
    end for;
    printf "SELFC %o %o %o\n", k, RealField(20)!Re(tot), RealField(20)!Im(tot);
end for;
R50 := RealField(50);
for key in Sort([k : k in Keys(Tg)]) do
    printf "TG %o g=%o %o %o\n", key[2], key[1], R50!Re(Tg[key]), R50!Im(Tg[key]);
end for;
printf "DONE\n";
quit;
