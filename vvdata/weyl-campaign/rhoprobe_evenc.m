// Probe the SHAPE of the even-c failure: rhokubota3 proved the fiber phases
// are coherent (sum = single term), so for the v2(c)=1 bad words the true
// rho(A)e_0 differs from the Stromberg single-term by something z-DEPENDENT.
// Dump, for chosen words, the ratio rv[z]/sv[z] over the support: magnitude,
// phase as a multiple of 1/levD, vs Qbar[z], the order of z, and the pairings
// B(z, w) with the 2-torsion elements w.  Also compare supports exactly.
//
//   magma -b DD:=15 NN:=2 vvdata/weyl-campaign/rhoprobe_evenc.m
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

// 2-torsion elements
tor2 := [ i : i in [1..n] | IsZero(2*elts[i]) ];
ordelt := function(z)
    k := 1; while not IsZero(k*elts[z]) do k +:= 1; end while; return k;
end function;
printf "2-torsion: %o elements\n", #tor2;

reps := VVCosetReps(M);
words := [ VVSTWord(g) : g in reps ];

DcTab := AssociativeArray();
XcTab := AssociativeArray();
getxc := procedure(cm, ~DcTab, ~XcTab)
    if IsDefined(DcTab, cm) then return; end if;
    dc := [ i : i in [1..n] | IsZero(cm*elts[i]) ];
    DcTab[cm] := dc;
    xcs := [ ];
    for i in [1..n] do
        ok := true;
        for j in dc do
            if frac(cm*Qbar[j] - BB(i, j)) ne 0 then ok := false; break; end if;
        end for;
        if ok then Append(~xcs, i); end if;
    end for;
    error if #xcs eq 0, "no x_c for", cm;
    XcTab[cm] := xcs;   // ALL solutions this time
end procedure;

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

R8 := RealField(8);
TOL := 10^(-30);

PROBE := [63, 126, 131, 3];   // three bad (c=-2,-6,-10) + one good control
for wi in PROBE do
    w := words[wi];
    gmat := VVWordMatrix(w);
    gi := gmat^(-1);
    a := gi[1][1]; b := gi[1][2]; c := gi[2][1]; d := gi[2][2];
    if c eq 0 then continue; end if;
    rv := VVRhoInvE0FFT(fftdata, w);
    cm := ((c mod levD) + levD) mod levD;
    getxc(cm, ~DcTab, ~XcTab);
    xcs := XcTab[cm]; dc := DcTab[cm];
    xc := xcs[1];
    printf "\n==== wi=%o (a,b,c,d)=(%o,%o,%o,%o)  #Dc=%o  #xc-solutions=%o\n",
        wi, a, b, c, d, #dc, #xcs;
    // how do the xc solutions distribute over cosets of c*D?
    // (two xc in the same coset of cD give the same formula up to global phase)
    cD := { idx[c*elts[j]] : j in [1..n] };
    printf "xc list (first 8): %o;  pairwise diffs in cD? %o\n",
        [ xcs[k] : k in [1..Minimum(8, #xcs)] ],
        [ idx[elts[xcs[k]] - elts[xc]] in cD : k in [2..Minimum(8, #xcs)] ];
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
    suppr := { i : i in [1..n] | Abs(rv[i]) gt TOL };
    supps := { i : i in [1..n] | Abs(sv[i]) gt TOL };
    printf "supports: |rv|=%o |sv|=%o  rv-only=%o sv-only=%o\n",
        #suppr, #supps, #(suppr diff supps), #(supps diff suppr);
    common := Sort([ z : z in suppr meet supps ]);
    // ratio structure over common support
    printf "z  60Q(z) ord(z)  |r|  levD*argfrac(r)  B(z,t2)*levD...\n";
    seen := AssociativeArray();
    for k in [1..#common] do
        z := common[k];
        r := rv[z]/sv[z];
        ph := frac(Arg(r)/(2*Pi(RealField(40))));
        key := <Integers()!(levD*Qbar[z]) mod levD,
                Round(levD*ph) mod levD>;
        if IsDefined(seen, key) then seen[key] +:= 1; else seen[key] := 1; end if;
        if k le 24 then
            printf "z=%o 60Q=%o ord=%o |r|=%o ph*%o=%o B2=%o\n",
                z, Integers()!(levD*Qbar[z]) mod levD, ordelt(z),
                R8!Abs(r), levD, R8!(levD*ph),
                [ Integers()!(levD*BB(z, t)) mod levD : t in tor2 ];
        end if;
    end for;
    printf "ratio census <60Q(z) mod %o, %o*argfrac(r)> -> count:\n", levD, levD;
    for k in Sort([ x : x in Keys(seen) ]) do
        printf "  %o -> %o\n", k, seen[k];
    end for;
end for;
quit;
