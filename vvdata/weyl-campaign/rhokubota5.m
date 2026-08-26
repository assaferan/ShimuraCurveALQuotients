// The even-c phase law, closed form: the rhoprobe2/rhokubota4 chain showed the
// even-c defect of the naive single-term formula is the character
// e((1-a)/2 B(t0, z)), and at 15_2 the x_c table satisfies x_c = -(c/2) t0
// for every v2(c) = 1 residue, which collapses the correction to putting an
// "a" on the B-term:  phase = e(a c Q(x) + a B(x_c, x))  for EVEN c (odd c
// keeps the plain B; both reduce to the same thing whenever (a-1) B(x_c, .)
// is integral, which is why 10_7 never showed the defect).  This script runs
// the aB variant and re-scores the unit u against the Kubota cocycle sign.
//
//   magma -b DD:=15 NN:=2 vvdata/weyl-campaign/rhokubota5.m
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

reps := VVCosetReps(M);
words := [ VVSTWord(g) : g in reps ];
printf "#cosets = %o\n", #words;
cvals := {};
for w in words do
    gmat := VVWordMatrix(w); gi := gmat^(-1);
    if gi[2][1] ne 0 then Include(~cvals, ((gi[2][1] mod levD) + levD) mod levD); end if;
end for;

DcTab := AssociativeArray();
XcTab := AssociativeArray();
for cm in cvals do
    dc := [ i : i in [1..n] | IsZero(cm*elts[i]) ];
    DcTab[cm] := dc;
    xc := 0;
    for i in [1..n] do
        ok := true;
        for j in dc do
            if frac(cm*Qbar[j] - BB(i, j)) ne 0 then ok := false; break; end if;
        end for;
        if ok then xc := i; break; end if;
    end for;
    error if xc eq 0, "no x_c for", cm;
    XcTab[cm] := xc;
end for;
printf "Dc/xc tables built (%o residues)\n", #cvals;

gsum := function(d, xi_)
    S := &+[ CC | ee(frac(d*Qbar[j] + (xi_ eq 0 select 0 else BB(xi_, j)))) : j in [1..n] ];
    nd := #[ i : i in [1..n] | IsZero(d*elts[i]) ];
    return S / Sqrt(CC!n * nd);
end function;
g3 := gsum(-1, 0)^3;

hilb := func< x, y | (x lt 0 and y lt 0) select -1 else 1 >;

epsDX := function(a, c)
    c2 := c; while IsEven(c2) do c2 := c2 div 2; end while;
    v := CC!KroneckerSymbol(-a, AbsoluteValue(c)) * hilb(c, -a);
    if IsEven(c) then v *:= ee(CC!((c2 + 1)*(a + 1))/8); end if;
    return v;
end function;

mn_of := function(a, b, c, d)
    for ntry := -4*levD to 4*levD do
        dp := c*ntry - d;
        if dp le 0 or GCD(dp, levD) ne 1 then continue; end if;
        for mtry := 1 to 8*levD do
            if (dp*mtry - c) mod levD ne 0 then continue; end if;
            if IsOdd(c) then
                if (dp - 1) mod 8 ne 0 then continue; end if;
                if (mtry - c) mod 8 ne 0 then continue; end if;
            else
                if (a*ntry - b) mod 8 ne 0 then continue; end if;
                if (mtry + a*c) mod 8 ne 0 then continue; end if;
            end if;
            return mtry, ntry;
        end for;
    end for;
    error "no (m, n) for", a, b, c, d;
end function;

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

R10 := RealField(10);
arg8 := func< z | R10!(8*Arg(z)/(2*Pi(RealField(30)))) >;
oddpart := function(c) c2 := AbsoluteValue(c); while IsEven(c2) do c2 := c2 div 2; end while; return c2; end function;

nhit := 0; nanti := 0; nbad := 0; nskip := 0;
printf "BAD-word table: wi (a,b,c,d) | v2(c) codd%%8 a%%8 d%%8 | arg8(u) |u| u8dev | propdev | ksign conv | arg8(u/ks)\n";
for wi->w in words do
    gmat := VVWordMatrix(w);
    gi := gmat^(-1);
    a := gi[1][1]; b := gi[1][2]; c := gi[2][1]; d := gi[2][2];
    if c eq 0 then nskip +:= 1; continue; end if;
    rv := VVRhoInvE0FFT(fftdata, w);
    cm := ((c mod levD) + levD) mod levD;
    xc := XcTab[cm]; dc := DcTab[cm];
    mm, nn := mn_of(a, b, c, d);
    dp := c*nn - d;
    xi := epsDX(a, c) * g3 * gsum(mm, xc eq i0 select 0 else xc) * gsum(dp, 0);
    fac := xi * Sqrt(CC!#dc / CC!n);
    aB := IsEven(c) select a else 1;   // the even-c law: a on the B-term
    sv := [ CC!0 : i in [1..n] ];
    for j in [1..n] do
        z := idx[c*elts[j] + elts[xc]];
        if sv[z] eq 0 then
            sv[z] := fac * ee(frac(a*c*Qbar[j] + (xc eq i0 select 0 else aB*BB(xc, j))));
        end if;
    end for;
    mx := 1;
    for i in [1..n] do if Abs(rv[i]) gt Abs(rv[mx]) then mx := i; end if; end for;
    if Abs(sv[mx]) lt 10^(-30) then
        printf "SUPPORT wi=%o (a,b,c,d)=(%o,%o,%o,%o)\n", wi, a, b, c, d;
        nskip +:= 1; continue;
    end if;
    u := rv[mx]/sv[mx];
    dev := Maximum([ Abs(rv[i] - u*sv[i]) : i in [1..n] ]);
    ksign, acc := kubota(w);
    error if acc ne gi and acc ne -gi, "word matrix mismatch", wi;
    conv := acc eq gi select 1 else -1;
    if dev lt 10^(-40) and Abs(u - ksign) lt 10^(-30) then
        nhit +:= 1;
    elif dev lt 10^(-40) and Abs(u + ksign) lt 10^(-30) then
        nanti +:= 1;
        printf "ANTI wi=%o (a,b,c,d)=(%o,%o,%o,%o) v2c=%o codd8=%o a8=%o d8=%o ks=%o conv=%o\n",
            wi, a, b, c, d, Valuation(c, 2), oddpart(c) mod 8,
            ((a mod 8)+8) mod 8, ((d mod 8)+8) mod 8, ksign, conv;
    else
        nbad +:= 1;
        printf "BAD wi=%o (%o,%o,%o,%o) | %o %o %o %o | %o %o %o | %o | %o %o | %o\n",
            wi, a, b, c, d,
            Valuation(c, 2), oddpart(c) mod 8, ((a mod 8)+8) mod 8, ((d mod 8)+8) mod 8,
            arg8(u), R10!Abs(u), R10!Abs(u^8 - 1),
            R10!dev,
            ksign, conv,
            arg8(u/ksign);
    end if;
end for;
printf "SUMMARY: HIT %o  ANTI %o  BAD %o  skipped %o (of %o)\n", nhit, nanti, nbad, nskip, #words;
quit;
