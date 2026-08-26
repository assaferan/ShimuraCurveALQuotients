// The tall-M oracle test: predict the banked 210_1 RHOV dump (960 words at
// M = 420, |D| = 88200) from the CLOSED FORM -- no FFT applies:
//   rv[z] = u(w) xi(a,c) sqrt(|D[c]|/|D|) e(ac Q(y) + B(x_c, y)),
//   z = c y + x_c,  u(w) = Kubota cocycle sign (bare-S word: extra -1),
//   x_c per Def 2.15 (Jordan-basis sum of the 2-part at v2(c) = 1, else 0),
//   xi(a,c) per Def 6.1 (per-Jordan-component Gauss sums).
// At 210_1 the dump tracks the ZERO coset (est = i0, #iso = 1), so per word
// we predict the coefficient of e_0: need y with c y = -x_c; if none, 0.
// Output: "PRED wi re im" on the same 50-digit grid; diff offline vs RHOV.
//
//   magma -b vvdata/weyl-campaign/rhokubota9_210.m
AttachSpec("ShimuraQuotients.spec");

D := 210; N := 1;

PREC := 80;
CC := ComplexField(PREC); ii := CC.1; pi := Pi(CC);
ee := func< z | Exp(2*pi*ii*z) >;

Ld := ShimuraCurveLattice(D, N);
M := IsOdd(D*N) select 4*D*N else 2*D*N;
fftdata := VVWeilFFT(Ld, CC : Dual := true);
elts := fftdata[7]; i0 := fftdata[8];
n := #elts;
printf "BASE %o %o  M = %o  |D| = %o  i0 = %o\n", D, N, M, n, i0;

Qr := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
lift := [ ChangeRing(elts[i]@@Ld`to_disc, Rationals()) : i in [1..n] ];
frac := func< r | r - Floor(r) >;
Qbar := [ frac(-(lift[i]*Qr, lift[i])/(2*dn^2)) : i in [1..n] ];
BB := function(i, j)
    return frac(-(lift[i]*Qr, lift[j])/dn^2);
end function;
idx := AssociativeArray();
for i in [1..n] do idx[elts[i]] := i; end for;

levD := M;
error if not forall{ j : j in [1..n] | frac(levD*Qbar[j]) eq 0 }, "levD != M";
printf "level(Dbar) = %o (verified)\n", levD;

// Sylow (Jordan-prime) index lists via p^3-annihilation (exponents here are
// squarefree-times-2, so p^3 is safely past every block)
plist := PrimeDivisors(n);
Syl := AssociativeArray();
for p in plist do
    Syl[p] := [ i : i in [1..n] | IsZero(p^3*elts[i]) ];
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

// Jordan basis of the elementary-abelian 2-part; canonical x_c (Def 2.15)
tor2nz := [ i : i in Syl[2] | not IsZero(elts[i]) and IsZero(2*elts[i]) ];
r2 := Valuation(#[ i : i in Syl[2] | IsZero(2*elts[i]) ], 2);
error if #tor2nz ne 2^r2 - 1, "2-torsion count mismatch";
error if r2 ne 3, "unexpected 2-rank";
basis2 := [ ]; found := false;
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
error if not found, "no Jordan basis for J_2";
xc2 := elts[basis2[1]] + elts[basis2[2]] + elts[basis2[3]];
xc2i := idx[xc2];
printf "Jordan basis J_2: %o  Q*4: %o  x_c = %o\n",
    basis2, [ Integers()!(4*Qbar[i]) mod 4 : i in basis2 ], xc2i;

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

// per-c-residue: a preimage y0 with c*y0 = -x_c (or 0 = no e_0 support), and |D[c]|
YTab := AssociativeArray();  // cm -> <y0 index or 0, dcn>
R50 := RealField(50);
for wi->w in words do
    gmat := VVWordMatrix(w);
    gi := gmat^(-1);
    a := gi[1][1]; b := gi[1][2]; c := gi[2][1]; d := gi[2][2];
    if c eq 0 then printf "PRED %o SKIP c0\n", wi; continue; end if;
    k2 := Valuation(c, 2);
    xci := (k2 eq 1) select xc2i else i0;
    cm := ((c mod levD) + levD) mod levD;
    key := <cm, xci>;
    if not IsDefined(YTab, key) then
        tgt := xci eq i0 select Parent(elts[1])!0 else -elts[xci];
        y0 := 0;
        dcn := 0;
        for j in [1..n] do
            if IsZero(cm*elts[j]) then dcn +:= 1; end if;
            if y0 eq 0 and cm*elts[j] eq tgt then y0 := j; end if;
        end for;
        YTab[key] := <y0, dcn>;
    end if;
    y0 := YTab[key][1]; dcn := YTab[key][2];
    if y0 eq 0 then
        printf "PRED %o ZERO\n", wi;
        continue;
    end if;
    // xi(a, c) per Def 6.1
    xiJ := CC!1;
    for p in plist do
        if c mod p ne 0 then
            xiJ *:= gsumJ(c, 0, p);
        else
            xiJ *:= CC!KroneckerSymbol(-a, #Syl[p]) *
                    gsumJ(-a*c, xci eq i0 select 0 else xci, p);
        end if;
    end for;
    xi0 := CC!KroneckerSymbol(-a, AbsoluteValue(c)) * hilb(-a, c);
    xi2 := IsOdd(c) select CC!1 else ee(CC!(-(a + 1)*(oddpart(c) - 1 + sign2))/8);
    e4s := ee(CC!(-signD)/4);
    xidef := e4s * xi0 * xi2 * xiJ;
    ksign, acc := kubota(w);
    error if acc ne gi and acc ne -gi, "word matrix mismatch", wi;
    u := ksign;
    if #w eq 1 and w[1][1] eq "S" then u := -u; end if;   // the bare-S rule
    pred := u * xidef * Sqrt(CC!dcn/CC!n) *
            ee(frac(a*c*Qbar[y0] + (xci eq i0 select 0 else BB(xci, y0))));
    printf "PRED %o %o %o\n", wi, R50!Re(pred), R50!Im(pred);
end for;
printf "DONE\n";
quit;
