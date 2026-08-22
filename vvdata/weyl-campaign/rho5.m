// Stromberg's explicit formula (Math. Z. 275 (2013), Thm 6.4) for the Weil
// representation of the module D-bar = (A, -q), verified coset by coset against
// the package's VVRhoInvE0FFT -- the first step of the theta_g phase-law proof.
// All caches precomputed at top level (Magma closures cannot mutate captures).
//
//   magma -b DD:=15 NN:=2 rho5.m
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

// level of Dbar
levD := 1;
for L in Divisors(2*M^2) do
    if forall{ j : j in [1..n] | frac(L*Qbar[j]) eq 0 } then levD := L; break; end if;
end for;
printf "level(Dbar) = %o\n", levD;

// coset words and the c-residues needed
reps := VVCosetReps(M);
words := [ VVSTWord(g) : g in reps ];
printf "#cosets = %o\n", #words;
cvals := {};
for w in words do
    gmat := VVWordMatrix(w); gi := gmat^(-1);
    if gi[2][1] ne 0 then Include(~cvals, ((gi[2][1] mod levD) + levD) mod levD); end if;
end for;

// precomputed D[c] and x_c tables keyed by c mod levD
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

// pure Gauss sum (no mutation)
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

// choose m, n per (6.1)-(6.3)
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

// ---- main loop ----
npass := 0; nfail := 0;
classsigns := AssociativeArray();
for wi->w in words do
    gmat := VVWordMatrix(w);
    gi := gmat^(-1);
    a := gi[1][1]; b := gi[1][2]; c := gi[2][1]; d := gi[2][2];
    if c eq 0 then continue; end if;
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
    mx := 1;
    for i in [1..n] do if Abs(rv[i]) gt Abs(rv[mx]) then mx := i; end if; end for;
    if Abs(sv[mx]) lt 10^(-30) then
        printf "SUPPORT MISMATCH coset %o (c=%o)\n", wi, gmat[2][1] mod M;
        nfail +:= 1; continue;
    end if;
    u := rv[mx]/sv[mx];
    dev := Maximum([ Abs(rv[i] - u*sv[i]) : i in [1..n] ]);
    okunit := Abs(u^8 - 1) lt 10^(-30);
    if dev lt 10^(-40) and okunit then
        npass +:= 1;
        cls := GCD(gmat[2][1] mod M, M);
        au := Round(8*Arg(u)/(2*Pi(RealField(30)))) mod 8;
        key := <cls, au>;
        if IsDefined(classsigns, key) then classsigns[key] +:= 1;
        else classsigns[key] := 1; end if;
    else
        nfail +:= 1;
        printf "FAIL coset %o (c=%o,d=%o): dev %o unit8 %o u=(%o,%o)\n", wi, c, d,
            RealField(6)!dev, okunit, RealField(8)!Re(u), RealField(8)!Im(u);
    end if;
end for;
printf "STROMBERG vs FFT: %o pass, %o fail (up to a global 8th-root unit per word)\n",
    npass, nfail;
printf "word-unit distribution <class, arg8(u)> -> count:\n";
for k in Sort([x : x in Keys(classsigns)]) do
    printf "  %o -> %o\n", k, classsigns[k];
end for;

// ---- canonical-rep breakdown per live class: every factor of xi, the word unit,
// and the measured eta*-entry ----
isoidx := [ i : i in [1..n] | Qbar[i] eq 0 ];
est := isoidx[2];
R10 := RealField(10);
arg8 := func< z | R10!(8*Arg(z)/(2*Pi(RealField(30)))) >;
for wi->w in words do
    gmat := VVWordMatrix(w);
    cw := gmat[2][1] mod M; dw := gmat[2][2] mod M;
    g := GCD(cw, M);
    if cw ne g or dw ne 1 then continue; end if;   // canonical rep (g, 1) only
    gi := gmat^(-1);
    a := gi[1][1]; b := gi[1][2]; c := gi[2][1]; d := gi[2][2];
    if c eq 0 then continue; end if;
    rv := VVRhoInvE0FFT(fftdata, w);
    if Abs(rv[est]) lt 10^(-40) then
        printf "CANON g=%o DEAD at eta*\n", g; continue;
    end if;
    cm := ((c mod levD) + levD) mod levD;
    xc := XcTab[cm]; dc := DcTab[cm];
    mm, nn := mn_of(a, b, c, d);
    dp := c*nn - d;
    epp := epsDX(a, c);
    gm := gsum(mm, xc eq i0 select 0 else xc);
    gd := gsum(dp, 0);
    xi := epp * g3 * gm * gd;
    // predicted eta*-entry: xi sqrt(|D[c]|/|D|) e(ac Qbar(y0) + B(xc,y0)), y0 fiber elt
    y0 := 0;
    for j in [1..n] do
        if idx[c*elts[j] + elts[xc]] eq est then y0 := j; break; end if;
    end for;
    pred := CC!0;
    if y0 ne 0 then
        pred := xi * Sqrt(CC!#dc/CC!n) *
                ee(frac(a*c*Qbar[y0] + (xc eq i0 select 0 else BB(xc, y0))));
    end if;
    u := Abs(pred) gt 10^(-40) select rv[est]/pred else CC!0;
    printf "CANON g=%o (a,b,c,d)=(%o,%o,%o,%o) m=%o n=%o dp=%o\n", g, a, b, c, d, mm, nn, dp;
    printf "  meas arg8 = %o   |meas| = %o\n", arg8(rv[est]), R10!Abs(rv[est]);
    printf "  u = (%o, %o) arg8 %o\n", R10!Re(u), R10!Im(u), arg8(u);
    printf "  eps = (%o, %o) arg8 %o;  g3 arg8 %o;  gm arg8 %o;  gd arg8 %o;  xi arg8 %o\n",
        R10!Re(epp), R10!Im(epp), arg8(epp), arg8(g3), arg8(gm), arg8(gd), arg8(xi);
    printf "  fiberphase arg8 %o\n",
        y0 eq 0 select R10!999 else
        arg8(ee(frac(a*c*Qbar[y0] + (xc eq i0 select 0 else BB(xc, y0)))));
end for;
quit;
