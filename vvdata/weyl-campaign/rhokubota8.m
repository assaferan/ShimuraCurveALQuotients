// Stromberg Thm 6.4 implemented VERBATIM, with the canonical x_c of
// Definition 2.15 (the piece every previous round missed):
//   x_c = 2^{k-1} sum(gamma_i) in the 2^k-Jordan component (k = ord_2(c)),
//   where {gamma_i} is a JORDAN (B-orthogonal) basis of that component;
//   x_c = 0 for odd c, and x_c = 0 when D has no 2^k-component -- which is
//   exactly why v2(c) >= 2 words never showed a defect at 15_2.
// xi(a,c) per Definition 6.1 (per-Jordan-component Gauss sums).  Column x=0:
//   rho(A) e_0 [cy + x_c] = xi(a,c) sqrt(|D[c]|/|D|) e(ac Q(y) + B(x_c, y)).
// Score u := rv/sv against the Kubota cocycle sign, with full-vector dev.
//
//   magma -b DD:=15 NN:=2 vvdata/weyl-campaign/rhokubota8.m
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

gsumJ := function(r, xi_, p)
    J := Syl[p];
    S := &+[ CC | ee(frac(r*Qbar[j] + (xi_ eq 0 select 0 else BB(xi_, j)))) : j in J ];
    nk := #[ j : j in J | IsZero(r*elts[j]) ];
    return S / Sqrt(CC!(#J) * nk);
end function;

R30 := RealField(30);
sig8 := function(p)
    g := gsumJ(1, 0, p);
    return Round(8*Arg(g)/(2*Pi(R30))) mod 8;
end function;
signD := (&+[ Integers() | sig8(p) : p in plist ]) mod 8;
sign2 := 2 in plist select sig8(2) else 0;
printf "sign(Dbar) = %o mod 8 (2-part %o)\n", signD, sign2;

// ---- Jordan basis of the 2-part of exponent-2 type (A_2 blocks) and the
// canonical x_c of Def 2.15 for k = 1.  (Bases here have (Z/2)^r 2-parts; a
// 2^k component with k >= 2 does not exist, so x_c = 0 for v2(c) >= 2.)
tor2nz := [ i : i in Syl[2] | ords[i] eq 2 ];
r2 := Valuation(#Syl[2], 2);
error if #tor2nz ne 2^r2 - 1, "2-part not elementary abelian; extend Def 2.15 handling";
// find a B-orthogonal basis by brute force among the nonzero 2-torsion
basis2 := [ ];
found := false;
if r2 eq 3 then
    for i1 in tor2nz do
        for i2 in tor2nz do
            if i2 eq i1 then continue; end if;
            if frac(BB(i1, i2)) ne 0 then continue; end if;
            for i3 in tor2nz do
                if i3 eq i1 or i3 eq i2 then continue; end if;
                if elts[i3] eq elts[i1] + elts[i2] then continue; end if;
                if frac(BB(i1, i3)) ne 0 or frac(BB(i2, i3)) ne 0 then continue; end if;
                basis2 := [i1, i2, i3]; found := true; break;
            end for;
            if found then break; end if;
        end for;
        if found then break; end if;
    end for;
else
    error "unexpected 2-rank", r2;
end if;
error if not found, "no Jordan basis found for the 2-part";
xc2 := elts[basis2[1]] + elts[basis2[2]] + elts[basis2[3]];
xc2i := idx[xc2];
printf "Jordan basis of J_2: %o  Q-values *4: %o\n",
    basis2, [ Integers()!(4*Qbar[i]) mod 4 : i in basis2 ];
printf "x_c (k=1) = index %o coords %o\n", xc2i, Eltseq(xc2);

reps := VVCosetReps(M);
words := [ VVSTWord(g) : g in reps ];
printf "#cosets = %o\n", #words;

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

nhit := 0; nanti := 0; nbad := 0; nskip := 0;
census := AssociativeArray();
for wi->w in words do
    gmat := VVWordMatrix(w);
    gi := gmat^(-1);
    a := gi[1][1]; b := gi[1][2]; c := gi[2][1]; d := gi[2][2];
    if c eq 0 then nskip +:= 1; continue; end if;
    rv := VVRhoInvE0FFT(fftdata, w);
    // canonical x_c (Def 2.15): k = v2(c); k=0 or no 2^k-component -> 0
    k2 := Valuation(c, 2);
    xc := (k2 eq 1) select xc2i else i0;   // i0 = identity index (x_c = 0)
    dcn := #[ i : i in [1..n] | IsZero(c*elts[i]) ];
    // xi(a, c) per Def 6.1
    xiJ := CC!1;
    for p in plist do
        if c mod p ne 0 then
            xiJ *:= gsumJ(c, 0, p);
        else
            xiJ *:= CC!KroneckerSymbol(-a, #Syl[p]) * gsumJ(-a*c, xc eq i0 select 0 else xc, p);
        end if;
    end for;
    xi0 := CC!KroneckerSymbol(-a, AbsoluteValue(c)) * hilb(-a, c);
    xi2 := IsOdd(c) select CC!1 else ee(CC!(-(a + 1)*(oddpart(c) - 1 + sign2))/8);
    e4s := ee(CC!(-signD)/4);
    xidef := e4s * xi0 * xi2 * xiJ;
    fac := xidef * Sqrt(CC!dcn / CC!n);
    sv := [ CC!0 : i in [1..n] ];
    for j in [1..n] do
        z := idx[c*elts[j] + (xc eq i0 select Parent(elts[1])!0 else elts[xc])];
        if sv[z] eq 0 then
            sv[z] := fac * ee(frac(a*c*Qbar[j] + (xc eq i0 select 0 else BB(xc, j))));
        end if;
    end for;
    mx := 1;
    for i in [1..n] do if Abs(rv[i]) gt Abs(rv[mx]) then mx := i; end if; end for;
    if Abs(sv[mx]) lt TOL then
        printf "SUPPORT wi=%o (a,b,c,d)=(%o,%o,%o,%o) k2=%o\n", wi, a, b, c, d, k2;
        nskip +:= 1; continue;
    end if;
    u := rv[mx]/sv[mx];
    dev := Maximum([ Abs(rv[i] - u*sv[i]) : i in [1..n] ]);
    ksign, acc := kubota(w);
    error if acc ne gi and acc ne -gi, "word matrix mismatch", wi;
    if dev lt 10^(-40) and Abs(u - ksign) lt TOL then
        nhit +:= 1;
    elif dev lt 10^(-40) and Abs(u + ksign) lt TOL then
        nanti +:= 1;
        printf "ANTI wi=%o (a,b,c,d)=(%o,%o,%o,%o) k2=%o ks=%o\n", wi, a, b, c, d, k2, ksign;
    else
        nbad +:= 1;
        printf "BAD wi=%o (%o,%o,%o,%o) k2=%o arg8(u)=%o |u|=%o dev=%o ks=%o\n",
            wi, a, b, c, d, k2, arg8(u), R8!Abs(u), R8!dev, ksign;
    end if;
    kk := <Round(60*Arg(u/ksign)/(2*Pi(R30))) mod 60, dev lt 10^(-40) select 1 else 0>;
    if IsDefined(census, kk) then census[kk] +:= 1; else census[kk] := 1; end if;
end for;
printf "SUMMARY: HIT %o  ANTI %o  BAD %o  skipped %o (of %o)\n", nhit, nanti, nbad, nskip, #words;
printf "census <60*arg(u/ks), propOK> -> count:\n";
for k in Sort([x : x in Keys(census)]) do printf "  %o -> %o\n", k, census[k]; end for;
quit;
