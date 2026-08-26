// Fit the even-c defect as a linear character: for every bad word at 15_2,
// find (w, u0) with  rv[z]/sv[z] = u0 * e(B(w, z))  on the support.  The
// x_c solution set is a full coset x_c + cD, and an in-coset shift t changes
// the single-term formula by e(c(1-a)B(t, .)), so the defect of an arbitrary
// x_c choice is exactly such a character; the fit should give w = (1-a)*t0(c)
// for a canonical t0 depending only on c (and the module).  u0 after the
// character is removed is the true word unit -- compare against Kubota.
//
//   magma -b DD:=15 NN:=2 vvdata/weyl-campaign/rhoprobe2_evenc.m
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

reps := VVCosetReps(M);
words := [ VVSTWord(g) : g in reps ];

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

R8 := RealField(8);
arg8 := func< z | R8!(8*Arg(z)/(2*Pi(RealField(30)))) >;
TOL := 10^(-30);

for wi->w in words do
    gmat := VVWordMatrix(w);
    gi := gmat^(-1);
    a := gi[1][1]; b := gi[1][2]; c := gi[2][1]; d := gi[2][2];
    if c eq 0 or IsOdd(c) then continue; end if;   // even-c words only
    rv := VVRhoInvE0FFT(fftdata, w);
    cm := ((c mod levD) + levD) mod levD;
    xc := XcTab[cm]; dc := DcTab[cm];
    mm, nn := mn_of(a, b, c, d);
    dp := c*nn - d;
    xi := epsDX(a, c) * g3 * gsum(mm, xc eq i0 select 0 else xc) * gsum(dp, 0);
    fac := xi * Sqrt(CC!#dc / CC!n);
    sv := [ CC!0 : i in [1..n] ];
    for j in [1..n] do
        z := idx[c*elts[j] + elts[xc]];
        if sv[z] eq 0 then
            sv[z] := fac * ee(frac(a*c*Qbar[j] + (xc eq i0 select 0 else BB(xc, j))));
        end if;
    end for;
    supp := [ i : i in [1..n] | Abs(rv[i]) gt TOL ];
    ksign, acc := kubota(w);
    // fit r(z) = u0 e(B(w0, z)); dedupe w0 by its character on the support
    found := [ ];
    for w0 in [1..n] do
        z1 := supp[1];
        u0 := (rv[z1]/sv[z1]) / ee(BB(w0, z1));
        ok := true;
        for k in [2..#supp] do
            z := supp[k];
            if Abs(rv[z]/sv[z] - u0*ee(BB(w0, z))) gt 10^(-25) then ok := false; break; end if;
        end for;
        if ok then
            Append(~found, <w0, u0>);
            if #found ge 3 then break; end if;   // a few witnesses is plenty
        end if;
    end for;
    if #found eq 0 then
        printf "wi=%o (a,c,d)=(%o,%o,%o) NO linear-character fit\n", wi, a, c, d;
    else
        w0 := found[1][1]; u0 := found[1][2];
        printf "wi=%o (a,c,d)=(%o,%o,%o) ks=%o  FIT w0=%o coords=%o arg8(u0)=%o u0^8dev=%o [%o wit]\n",
            wi, a, c, d, ksign, w0, Eltseq(elts[w0]), arg8(u0), R8!Abs(u0^8-1), #found;
    end if;
end for;
quit;
