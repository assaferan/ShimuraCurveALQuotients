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

// ---- kubota-cocycle prediction of the per-word unit u ----
// letters applied as inverses, g_1 first: operator = rho(g_k^-1)...rho(g_1^-1)
// = (accumulated cocycle) * rho(gamma_w^-1).  Kubota: sigma(g)sigma(h) =
// c(g,h) sigma(gh), c(g,h) = hilb(s(g)s(gh), s(h)s(gh)), s = (c ne 0 -> c, else d).
sfun := func< g | g[2][1] ne 0 select g[2][1] else g[2][2] >;
Smat := Matrix(Integers(), 2, 2, [0, -1, 1, 0]);
Tmat := func< k | Matrix(Integers(), 2, 2, [1, k, 0, 1]) >;
kubota := function(w)
    // matrices of the applied letters, in application order
    ls := [ ];
    for t in w do
        Append(~ls, t[1] eq "S" select Smat^(-1) else Tmat(-t[2]));
    end for;
    // note: whether T-letter inverts (Tdiag^(-k)) matches Tmat(-k) here
    acc := ls[1]; sign := 1;
    for i in [2..#ls] do
        h := ls[i]; prod := h*acc;   // operator composes on the LEFT
        sign *:= (sfun(h)*sfun(prod) lt 0 and sfun(acc)*sfun(prod) lt 0) select -1 else 1;
        acc := prod;
    end for;
    return sign, acc;
end function;

// ---- main loop: measure u via FFT (small base!), compare with kubota ----
npass := 0; nfail := 0; nskip := 0;
mism := AssociativeArray();
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
    sv := [ CC!0 : i in [1..n] ];
    for j in [1..n] do
        z := idx[c*elts[j] + elts[xc]];
        if sv[z] eq 0 then
            sv[z] := fac * ee(frac(a*c*Qbar[j] + (xc eq i0 select 0 else BB(xc, j))));
        end if;
    end for;
    mx := 1;
    for i in [1..n] do if Abs(rv[i]) gt Abs(rv[mx]) then mx := i; end if; end for;
    if Abs(sv[mx]) lt 10^(-30) then nskip +:= 1; continue; end if;
    u := rv[mx]/sv[mx];
    ksign, acc := kubota(w);
    error if acc ne gi and acc ne -gi, "word matrix mismatch", wi;
    // hypothesis: u = ksign * (convention factor depending on acc vs gi sign)
    conv := acc eq gi select 1 else -1;
    for hyp in [1..4] do
        pred := [ ksign, ksign*conv, conv, 1 ][hyp];
        key := <hyp, Abs(u - pred) lt 10^(-30) select "HIT" else (Abs(u + pred) lt 10^(-30) select "ANTI" else "OTHER")>;
        if IsDefined(mism, key) then mism[key] +:= 1; else mism[key] := 1; end if;
    end for;
    npass +:= 1;
    if wi mod 16 eq 0 then printf "...word %o\n", wi; end if;
end for;
printf "words tested %o (skipped %o)\n", npass, nskip;
printf "hypothesis scorecard (hyp, verdict) -> count:\n";
for k in Sort([x : x in Keys(mism)]) do printf "  %o -> %o\n", k, mism[k]; end for;
quit;
