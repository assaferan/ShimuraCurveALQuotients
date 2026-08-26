// Stromberg Thm 6.4 / Definition 6.1 DIRECT test.  At x = 0 the theorem reads
//   rho(A) e_0 [c y + x_c] = xi(a,c) sqrt(|D[c]|/|D|) e(ac Q(y) + B(x_c, y)),
// with xi(a,c) = e4(-sign(D)) xi0 xi2 prod_J xi(J) over nontrivial Jordan
// components (Def 6.1):
//   xi(J) = G(c, 0; J)                       if p not| c,
//   xi(J) = (-a/|J|) G(-ac, x_c; J)          if p | c,
//   xi0 = (-a/c)(-a,c)_oo, xi2 = 1 (c odd) or e8(-(a+1)(c2-1+sign(D_2)))
// (odd signature; even signature xi0 = xi2 = 1).
// The FFT operator differs from rho(A) by the empirically-pinned aB shape
// twist for even c; so measure  xi_meas := rv[z]/(sqrt(|D[c]|/|D|) *
// e(acQ(j) + aB B(x_c,j)))  (exact, shape-verified) and dump the ratio
// xi_meas / xi_def61 per word.  If the ratio is a pure 8th root for ALL
// words (incl. the 17 residual ones), Def 6.1's per-component Gauss sums
// explain the odd-part constants and only convention bookkeeping remains.
//
//   magma -b DD:=15 NN:=2 vvdata/weyl-campaign/rhokubota6.m
AttachSpec("ShimuraQuotients.spec");

D := 15; N := 2;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;

PREC := 80;
CC := ComplexField(PREC); ii := CC.1; pi := Pi(CC);
ee := func< z | Exp(2*pi*ii*z) >;

Ld := ShimuraCurveLattice(D, N);
M := IsOdd(D*N) select 4*D*N else 2*D*N;
fftdata := VVWeilFFT(Ld, CC : Dual := true);
elts := fftdata[7]; i0 := fftdata[8];
n := #elts;
printf "BASE %o %o  M = %o  |D| = %o\n", D, N, M, n;

Qr := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
lift := [ ChangeRing(elts[i]@@Ld`to_disc, Rationals()) : i in [1..n] ];
frac := func< r | r - Floor(r) >;
Qbar := [ frac(-(lift[i]*Qr, lift[i])/(2*dn^2)) : i in [1..n] ];
BB := function(i, j)
    return frac(-(lift[i]*Qr, lift[j])/dn^2);
end function;
idx := AssociativeArray();
for i in [1..n] do idx[elts[i]] := i; end for;

levD := 1;
for L in Divisors(2*M^2) do
    if forall{ j : j in [1..n] | frac(L*Qbar[j]) eq 0 } then levD := L; break; end if;
end for;
printf "level(Dbar) = %o\n", levD;

// element orders and p-Sylow (Jordan component) index lists
ordelt := function(z)
    k := 1; while not IsZero(k*elts[z]) do k +:= 1; end while; return k;
end function;
ords := [ ordelt(i) : i in [1..n] ];
plist := PrimeDivisors(n);
Syl := AssociativeArray();
for p in plist do
    Syl[p] := [ i : i in [1..n] | Set(PrimeDivisors(ords[i])) subset {p} ];
    printf "p = %o: |J_p| = %o\n", p, #Syl[p];
end for;

// normalized component Gauss sum G(r, x; J_p)
gsumJ := function(r, xi_, p)
    J := Syl[p];
    S := &+[ CC | ee(frac(r*Qbar[j] + (xi_ eq 0 select 0 else BB(xi_, j)))) : j in J ];
    nk := #[ j : j in J | IsZero(r*elts[j]) ];
    return S / Sqrt(CC!(#J) * nk);
end function;

// signatures: G(1, 0; J) = e8(sign(J))
R30 := RealField(30);
sig8 := function(p)
    g := gsumJ(1, 0, p);
    return Round(8*Arg(g)/(2*Pi(R30))) mod 8;
end function;
signD := (&+[ Integers() | sig8(p) : p in plist ]) mod 8;
sign2 := 2 in plist select sig8(2) else 0;
printf "sign(Dbar) = %o mod 8 (2-part %o)\n", signD, sign2;

reps := VVCosetReps(M);
words := [ VVSTWord(g) : g in reps ];
printf "#cosets = %o\n", #words;

DcTab := AssociativeArray();
XcTab := AssociativeArray();
for w in words do
    gmat := VVWordMatrix(w); gi := gmat^(-1);
    c := gi[2][1];
    if c eq 0 then continue; end if;
    cm := ((c mod levD) + levD) mod levD;
    if IsDefined(DcTab, cm) then continue; end if;
    dc := [ i : i in [1..n] | IsZero(cm*elts[i]) ];
    DcTab[cm] := dc;
    for i in [1..n] do
        ok := true;
        for j in dc do
            if frac(cm*Qbar[j] - BB(i, j)) ne 0 then ok := false; break; end if;
        end for;
        if ok then XcTab[cm] := i; break; end if;
    end for;
end for;

hilb := func< x, y | (x lt 0 and y lt 0) select -1 else 1 >;
oddpart := function(c) c2 := AbsoluteValue(c); while IsEven(c2) do c2 := c2 div 2; end while; return c2; end function;

sfun := func< g | g[2][1] ne 0 select g[2][1] else g[2][2] >;
Smat := Matrix(Integers(), 2, 2, [0, -1, 1, 0]);
Tmat := func< k | Matrix(Integers(), 2, 2, [1, k, 0, 1]) >;
kubota := function(w)
    ls := [ ];
    for t in w do
        Append(~ls, t[1] eq "S" select Smat^(-1) else Tmat(-t[2]));
    end for;
    acc := ls[1]; sign := 1;
    for i in [2..#ls] do
        h := ls[i]; prod := h*acc;
        sign *:= (sfun(h)*sfun(prod) lt 0 and sfun(acc)*sfun(prod) lt 0) select -1 else 1;
        acc := prod;
    end for;
    return sign, acc;
end function;

R8 := RealField(8);
arg8 := func< z | R8!(8*Arg(z)/(2*Pi(RealField(30)))) >;
TOL := 10^(-30);

census := AssociativeArray();
for wi->w in words do
    gmat := VVWordMatrix(w);
    gi := gmat^(-1);
    a := gi[1][1]; b := gi[1][2]; c := gi[2][1]; d := gi[2][2];
    if c eq 0 then continue; end if;
    rv := VVRhoInvE0FFT(fftdata, w);
    cm := ((c mod levD) + levD) mod levD;
    xc := XcTab[cm]; dc := DcTab[cm];
    aB := IsEven(c) select a else 1;
    // xi_meas from the max entry via the (shape-exact) aB formula
    mx := 1;
    for i in [1..n] do if Abs(rv[i]) gt Abs(rv[mx]) then mx := i; end if; end for;
    // find a preimage j with c*elts[j] + elts[xc] = elts[mx]
    jm := 0;
    for j in [1..n] do
        if idx[c*elts[j] + elts[xc]] eq mx then jm := j; break; end if;
    end for;
    error if jm eq 0, "no preimage", wi;
    shape := Sqrt(CC!#dc/CC!n) *
             ee(frac(a*c*Qbar[jm] + (xc eq i0 select 0 else aB*BB(xc, jm))));
    xim := rv[mx]/shape;
    // Def 6.1 product over components
    xiJ := CC!1;
    for p in plist do
        if c mod p ne 0 then
            xiJ *:= gsumJ(c, 0, p);
        else
            xiJ *:= CC!KroneckerSymbol(-a, #Syl[p]) * gsumJ(-a*c, xc, p);
        end if;
    end for;
    // odd-signature prefactors (Def 6.1)
    xi0 := CC!KroneckerSymbol(-a, AbsoluteValue(c)) * hilb(-a, c);
    xi2 := IsOdd(c) select CC!1 else ee(CC!(-(a + 1)*(oddpart(c) - 1 + sign2))/8);
    e4s := ee(CC!(-signD)/4);
    xidef := e4s * xi0 * xi2 * xiJ;
    ksign, acc := kubota(w);
    rat := xim/xidef;
    key := <arg8(rat), ksign>;
    kk := <Round(60*Arg(rat)/(2*Pi(R30))) mod 60, ksign, IsEven(c) select 1 else 0>;
    if IsDefined(census, kk) then census[kk] +:= 1; else census[kk] := 1; end if;
    if Abs(Abs(rat) - 1) gt 10^(-25) then
        printf "NONUNIT wi=%o (a,c,d)=(%o,%o,%o) |rat|=%o\n", wi, a, c, d, R8!Abs(rat);
    end if;
end for;
printf "census <60*arg(xi_meas/xi_def61), ksign, evenc> -> count:\n";
for k in Sort([x : x in Keys(census)]) do printf "  %o -> %o\n", k, census[k]; end for;
quit;
