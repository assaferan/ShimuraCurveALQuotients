// eis32s: eis32k with the O(|A|) setup replaced by a COORDINATE MODEL of the
// discriminant group.  Mathematically identical -- same Kubota closed form,
// same xi, same x_c -- but nothing ever enumerates A.
//
// A = (+)_r Z/m_r on the Smith generators of Ld`disc_grp, so an element is a
// k-tuple of integers (k <= 3).  Everything the driver needs is then O(1) or
// O(|J_p|) instead of O(|A|):
//   Qbar, B      -- from the k x k Gram on the generators (quadratic/bilinear)
//   |A[c]|       -- prod_r gcd(c, m_r)          (was an |A| scan)
//   c*y = tgt    -- k linear congruences        (was an |A| search)
//   J_p          -- built from the p-parts of the m_r, |J_p| <= 121 at 1155
//   the a=0 word -- closed form, see below      (was a full FFT over |A|)
// Dropped entirely: elts, lift, Qbar[], idx, the per-prime |A| Sylow scans.
//
// At D = 1155 this takes the setup from 2,668,050 objects to ~200.
//
//   magma -b DD:=1155 NN:=1 PR:=120 EF:=<pool> eis32s.m
// CTL:=1    word-by-word FFT control (small |A| only) -- the validation gate.
// VERIFY:=1 checks the coordinate model against the library's own element
//           list, elementwise (small |A| only).
AttachSpec("ShimuraQuotients.spec");

D := 15; N := 2;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;

PREC := 80;
if assigned PR then PREC := StringToInteger(PR); end if;
CC := ComplexField(PREC); ii := CC.1; pi := Pi(CC);
ee := func< z | Exp(2*pi*ii*z) >;

M := IsOdd(D*N) select 4*D*N else 2*D*N;
ds := Divisors(M); nd := #ds;
printf "BASE %o %o  M = %o  ds = %o\n", D, N, M, ds;
Ld := ShimuraCurveLattice(D, N);

// ---- the coordinate model ------------------------------------------------
G := Ld`disc_grp;
mods := Moduli(G);
keep := [ r : r in [1..#mods] | mods[r] gt 1 ];
ms := [ mods[r] : r in keep ];
k := #ms;
n := &*ms;
printf "|D| = %o  invariants = %o\n", n, ms;

Qr := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
frac := func< r | r - Floor(r) >;

// L = dn*Z^3 (RSpaceWithBasis(ScalarMatrix(3,denom))).  B(L, L^dual) is
// automatically integral, so Qbar and B descend to A -- hence are computable
// from the generators alone -- exactly when Q has even diagonal.  Check it:
// without this the coordinate model would be ill-defined.
for aa in [1..3] do
    error if IsOdd(Integers()!Qr[aa][aa]),
        "Q has odd diagonal: Qbar does not descend, coordinate model invalid", aa;
end for;

gens := [ G.(keep[r]) : r in [1..k] ];
wg := [ ChangeRing(g@@Ld`to_disc, Rationals()) : g in gens ];
QG := [ frac(-(wg[r]*Qr, wg[r])/(2*dn^2)) : r in [1..k] ];
BG := [ [ frac(-(wg[r]*Qr, wg[s])/dn^2) : s in [1..k] ] : r in [1..k] ];

zero := [ Integers() | 0 : r in [1..k] ];
addc := func< c, e | [ (c[r] + e[r]) mod ms[r] : r in [1..k] ] >;
subc := func< c, e | [ (c[r] - e[r]) mod ms[r] : r in [1..k] ] >;
mulc := func< t, c | [ (t*c[r]) mod ms[r] : r in [1..k] ] >;
iszc := func< c | forall{ r : r in [1..k] | c[r] mod ms[r] eq 0 } >;

Qc := function(c)
    s := &+[ Rationals() | c[r]^2*QG[r] : r in [1..k] ];
    for r in [1..k] do
        for s2 in [r+1..k] do s +:= c[r]*c[s2]*BG[r][s2]; end for;
    end for;
    return frac(s);
end function;
Bc := function(c, e)
    s := Rationals()!0;
    for r in [1..k] do
        for s2 in [1..k] do s +:= c[r]*e[s2]*BG[r][s2]; end for;
    end for;
    return frac(s);
end function;

levD := M;
error if exists{ r : r in [1..k] | frac(levD*QG[r]) ne 0 }, "levD != M (Q part)";
error if exists{ <r,s2> : r, s2 in [1..k] | frac(levD*BG[r][s2]) ne 0 }, "levD != M (B part)";
printf "level(Dbar) = %o (verified on the generators)\n", levD;

// ---- Sylow (Jordan-prime) components, built directly ---------------------
plist := PrimeDivisors(n);
JP := AssociativeArray();   // p -> list of coordinate vectors of J_p
LEVP := AssociativeArray(); // p -> modulus in which the multiplier r matters
for p in plist do
    cur := [ zero ];
    for r in [1..k] do
        ar := Valuation(ms[r], p);
        step := ms[r] div p^ar;
        nxt := [ ];
        for c in cur do
            for t in [0..p^ar - 1] do
                cc := c; cc[r] := (t*step) mod ms[r]; Append(~nxt, cc);
            end for;
        end for;
        cur := nxt;
    end for;
    JP[p] := cur;
    // r enters only through r*Qbar(j) mod 1 and r*j, so reducing r modulo the
    // lcm of the exponent and the Qbar denominators is exact.  (At p = 2 the
    // form has denominator 4 on a Z/2, so the exponent alone is NOT enough.)
    ex := &*[ Integers() | p^Valuation(ms[r], p) : r in [1..k] ];
    lp := LCM([ Integers() | ex ] cat [ Denominator(Qc(j)) : j in cur ]);
    LEVP[p] := lp;
    printf "p = %o: |J_p| = %o  (multiplier modulus %o)\n", p, #cur, lp;
end for;

// ---- Jordan basis of the elementary-abelian 2-part; canonical x_c --------
J2 := JP[2];
tor2 := [ j : j in J2 | iszc(mulc(2, j)) ];
tor2nz := [ j : j in tor2 | not iszc(j) ];
r2 := Valuation(#tor2, 2);
error if #tor2nz ne 2^r2 - 1, "2-torsion count mismatch";
error if r2 notin {1, 3}, "unexpected 2-rank", r2;
basis2 := [ ]; found := false;
if r2 eq 1 then
    // ODD DN: level 4DN, the 2-part is a single Z/2 of ODD type, so the
    // Jordan basis is its one generator and x_c is forced -- no search, no
    // ordering dependence.
    basis2 := [ tor2nz[1] ]; found := true;
    error if IsEven(Integers()!(4*Qc(tor2nz[1]))),
        "2-part is EVEN type: x_c rule not established here";
end if;
for i1 in (found select [ ] else tor2nz) do
    for i2 in tor2nz do
        if i2 eq i1 then continue; end if;
        if Bc(i1, i2) ne 0 then continue; end if;
        for i3 in tor2nz do
            if i3 eq i1 or i3 eq i2 then continue; end if;
            if i3 eq addc(i1, i2) then continue; end if;
            if Bc(i1, i3) ne 0 or Bc(i2, i3) ne 0 then continue; end if;
            basis2 := [i1, i2, i3]; found := true; break;
        end for;
        if found then break; end if;
    end for;
    if found then break; end if;
end for;
error if not found, "no Jordan basis for J_2";
xc := zero;
for b in basis2 do xc := addc(xc, b); end for;
printf "Jordan basis J_2 (2-rank %o): %o  Q*4: %o  x_c = %o\n",
    r2, basis2, [ Integers()!(4*Qc(b)) mod 4 : b in basis2 ], xc;

// ---- Gauss sums over the components --------------------------------------
// The multiplier enters only as r mod LEVP[p], and LEVP[p] <= 121 here, so
// tabulate every residue up front rather than caching inside a closure: the
// whole table is ~35k exponentials at 1155 and the word loop then reads it.
GTAB := AssociativeArray();   // p -> [ <g(r, xc=no), g(r, xc=yes)> : r ]
for p in plist do
    J := JP[p];
    QJp := [ Qc(j) : j in J ];
    BXCp := [ Bc(xc, j) : j in J ];
    row := [ ];
    for rr in [0..LEVP[p] - 1] do
        S0 := CC!0; S1 := CC!0; nk := 0;
        for jx in [1..#J] do
            t := rr*QJp[jx];
            S0 +:= ee(frac(t));
            S1 +:= ee(frac(t + BXCp[jx]));
            if iszc(mulc(rr, J[jx])) then nk +:= 1; end if;
        end for;
        nrm := Sqrt(CC!(#J) * nk);
        Append(~row, <S0/nrm, S1/nrm>);
    end for;
    GTAB[p] := row;
end for;
gsumJ := func< r, useXc, p |
    useXc select GTAB[p][(r mod LEVP[p]) + 1][2] else GTAB[p][(r mod LEVP[p]) + 1][1] >;

R30 := RealField(30);
sig8 := function(p)
    g := gsumJ(1, false, p);
    return Round(8*Arg(g)/(2*Pi(R30))) mod 8;
end function;
signD := (&+[ Integers() | sig8(p) : p in plist ]) mod 8;
sign2 := 2 in plist select sig8(2) else 0;
printf "sign(Dbar) = %o mod 8 (2-part %o)\n", signD, sign2;

// ---- the tracked coset ---------------------------------------------------
// At N = 1 there is exactly one isotropic coset (#iso = 2N-1), the zero one,
// so est is known without touching the group.  At N > 1 the convention is
// "the second isotropic coset in the library's element order", which is an
// ordering convention -- reproduce it by listing G (small bases only).
NEEDFULL := (N gt 1) or assigned CTL or assigned VERIFY
            or (assigned EST and StringToInteger(EST) ne 0);
estc := zero; esti := 0; elts := [ ]; coordl := [ ]; i0 := 0;
if NEEDFULL then
    elts := [ g : g in G ];
    error if #elts ne n, "element enumeration disagrees with |A|";
    coordl := [ [ (Eltseq(e)[keep[r]]) mod ms[r] : r in [1..k] ] : e in elts ];
    i0 := rep{ i : i in [1..n] | IsZero(elts[i]) };
    isoidx := [ i : i in [1..n] | Qc(coordl[i]) eq 0 ];
    esti := #isoidx ge 2 select isoidx[2] else isoidx[1];
    if assigned EST then
        estk := StringToInteger(EST);
        esti := estk eq 0 select i0 else isoidx[estk];
        printf "EST OVERRIDE: %o -> tracking coset index %o (i0 = %o)\n", estk, esti, i0;
    end if;
    estc := coordl[esti];
    printf "#isotropic = %o, tracking coset index %o  coords %o\n", #isoidx, esti, estc;
else
    printf "#isotropic = 1 (N = 1), tracking the ZERO coset -- group never enumerated\n";
end if;
// the a = 0 closed form below needs this
error if Qc(estc) ne 0, "tracked coset is not isotropic", estc;

// ---- VERIFY: the coordinate model against the library's own element list --
if assigned VERIFY then
    Qref := [ frac(-(ChangeRing(elts[i]@@Ld`to_disc, Rationals())*Qr,
                     ChangeRing(elts[i]@@Ld`to_disc, Rationals()))/(2*dn^2)) : i in [1..n] ];
    bad := 0;
    for i in [1..n] do
        if Qc(coordl[i]) ne Qref[i] then bad +:= 1; end if;
    end for;
    printf "VERIFY Qbar: %o/%o mismatches\n", bad, n;
    badB := 0;
    for i in [1..Min(n, 200)] do
        li := ChangeRing(elts[i]@@Ld`to_disc, Rationals());
        for j in [1..Min(n, 200)] do
            lj := ChangeRing(elts[j]@@Ld`to_disc, Rationals());
            if Bc(coordl[i], coordl[j]) ne frac(-(li*Qr, lj)/dn^2) then badB +:= 1; end if;
        end for;
    end for;
    printf "VERIFY B: %o mismatches on the 200x200 corner\n", badB;
    for p in plist do
        ref := #[ i : i in [1..n] | IsZero(p^3*elts[i]) ];
        printf "VERIFY J_%o: model %o  library %o  %o\n", p, #JP[p], ref,
            (#JP[p] eq ref) select "OK" else "MISMATCH";
    end for;
    error if bad ne 0 or badB ne 0, "coordinate model disagrees with the library";
end if;

reps := VVCosetReps(M);
words := [ VVSTWord(g) : g in reps ];
printf "#cosets = %o\n", #words;

hilb := func< x, y | (x lt 0 and y lt 0) select -1 else 1 >;
oddpart := function(c) c2 := AbsoluteValue(c); while IsEven(c2) do c2 := c2 div 2; end while; return c2; end function;

sfun := func< g | g[2][1] ne 0 select g[2][1] else g[2][2] >;
Smat := Matrix(Integers(), 2, 2, [0, -1, 1, 0]);
Tmat := func< kk | Matrix(Integers(), 2, 2, [1, kk, 0, 1]) >;
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

// ---- rho vector via the closed form --------------------------------------
nw := #words;
rhov := [ CC | 0 : wi in [1..nw] ];
// per-c-residue: a preimage y with c*y = est - x_c (or none), and |A[c]|
YTab := AssociativeArray();  // <cm, useXc> -> <ok, y, dcn>
for wi->w in words do
    if #w eq 0 or VVWordMatrix(w)[2][1] eq 0 then
        // T^k coset: rho(T^k)^{-1} e_0 = e_0, so the tracked component is 1 at
        // the ZERO coset and 0 at any nonzero one.
        rhov[wi] := iszc(estc) select CC!1 else CC!0;
        continue;
    end if;
    gmat := VVWordMatrix(w);
    gi := gmat^(-1);
    a := gi[1][1]; b := gi[1][2]; c := gi[2][1]; d := gi[2][2];
    if a eq 0 then
        // The S-coset, in closed form -- eis32k measured this one entry with a
        // full FFT over |A|, which is the last O(|A|) step in the driver.
        // Walk VVRhoInvE0FFT by hand on a word with a single S:
        //   * leading T^k act on e_0 by Tdiag[i0]^-k = e(-Q(0))^-k = 1;
        //   * the S sends e_0 to the constant vector conj(c_S) = e(1/8)/sqrt|A|
        //     (the transform of a delta at the origin is all-ones);
        //   * trailing T^k scale entry est by Tdiag[est]^-k = e(-Q(est))^-k = 1,
        //     since the tracked coset is ISOTROPIC.
        // So the entry is e(1/8)/sqrt|A| regardless of the T-letters and
        // regardless of which isotropic coset is tracked.  Matches the FFT
        // measurement on all eleven banked bases (|A| from 200 to 217800, both
        // parities of DN, both 2-ranks) -- including 15_2, where est is
        // NONZERO, which is what pins the coset-independence.
        error if #[ t : t in w | t[1] eq "S" ] ne 1,
            "a = 0 word has more than one S: closed form not established", wi;
        rhov[wi] := ee(CC!1/8) / Sqrt(CC!n);
        continue;
    end if;
    k2 := Valuation(c, 2);
    useXc := (k2 eq 1);
    cm := ((c mod levD) + levD) mod levD;
    key := <cm, useXc>;
    if not IsDefined(YTab, key) then
        tgt := useXc select subc(estc, xc) else estc;
        okY := true; y := [ Integers() | ];
        for r in [1..k] do
            m := ms[r]; g := GCD(cm, m);
            if tgt[r] mod g ne 0 then okY := false; break; end if;
            mm := m div g;
            if mm eq 1 then Append(~y, 0);
            else Append(~y, (InverseMod((cm div g) mod mm, mm) * ((tgt[r] div g) mod mm)) mod mm);
            end if;
        end for;
        dcn := &*[ Integers() | GCD(cm, ms[r]) : r in [1..k] ];
        YTab[key] := <okY, y, dcn>;
    end if;
    okY := YTab[key][1]; y := YTab[key][2]; dcn := YTab[key][3];
    if not okY then continue; end if;   // rhov[wi] stays 0
    // xi(a, c) per Def 6.1
    xiJ := CC!1;
    for p in plist do
        if c mod p ne 0 then
            xiJ *:= gsumJ(c, false, p);
        else
            xiJ *:= CC!KroneckerSymbol(-a, #JP[p]) * gsumJ(-a*c, useXc, p);
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
    rhov[wi] := u * xidef * Sqrt(CC!dcn/CC!n) *
            ee(frac(a*c*Qc(y) + (useXc select Bc(xc, y) else 0)));
    if wi mod 512 eq 0 then printf "rho %o/%o (closed form)\n", wi, nw; end if;
end for;
printf "rhov done (coordinate model, %o YTab keys)\n", #Keys(YTab);

// ---- CTL: word-by-word FFT control (affordable only at small |A|) --------
if assigned CTL then
    fftdata := VVWeilFFT(Ld, CC : Dual := true);
    error if fftdata[8] ne i0, "FFT zero index disagrees with the model";
    worst := RealField(30)!0; wworst := 0; nnz := 0;
    for wi->w in words do
        rvF := VVRhoInvE0FFT(fftdata, w);
        dev := Abs(rvF[esti] - rhov[wi]);
        if Abs(rvF[esti]) gt 10^(-30) then nnz +:= 1; end if;
        if dev gt worst then worst := RealField(30)!dev; wworst := wi; end if;
    end for;
    printf "CTL: worst |closed - FFT| = %o at word %o  (%o/%o FFT entries nonzero)\n",
        worst, wworst, nnz, nw;
    printf "CTL VERDICT: %o\n", worst lt 10^(-25) select "PASS" else "FAIL";
end if;

// ---- E pool (external only) ---------------------------------------------
cs := Divisors(M);
ordvec := function(r)  // 24*ord at cusp c, up to positive width factors
    return [ &+[ Rationals() | r[i]*GCD(c, ds[i])^2/ds[i] : i in [1..nd] ] : c in cs ];
end function;

Epool := {@ @};
if assigned EF then
    EXTRAE := eval Read(EF);
    for r in EXTRAE do
        rr := [ Integers() | x : x in r ];
        error if #rr ne nd, "EF vector length mismatch";
        error if exists{ o : o in ordvec(rr) | o lt 0 }, "EF vector not holomorphic";
        Include(~Epool, rr);
    end for;
    printf "#E's after external pool = %o\n", #Epool;
end if;
for Ei->rE in Epool do
    printf "EPOOL %o %o\n", Ei, rE;
end for;

// ---- cusp6 coset machinery, weight-general ------------------------------
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

SS := PowerSeriesRing(CC); t := SS.1;

// a0 of the eta monomial with exponent vector r (weight (sum r)/2) at coset w
a0at := function(w, r, wt)
    g := VVWordMatrix(w);
    tri := [ ];
    for d in ds do
        a, b, e, gd := triang(g, d);
        Append(~tri, <a, b, e, gd>);
    end for;
    W := LCM([ tri[i][3] : i in [1..nd] ]);
    L := &+[ Integers() | r[i]*tri[i][1]*(W div tri[i][3]) : i in [1..nd] ];
    if L gt 0 then return CC!0, L, W; end if;
    depth := -L + 1;
    units := [ ];
    for i->d in ds do
        a, b, e := Explode([tri[i][1], tri[i][2], tri[i][3]]);
        step := 24*a*(W div e);
        u := SS!1 + O(t^depth);
        nn := 1;
        while nn*step lt depth do
            u *:= 1 - ee(CC!(nn*b/e))*t^(nn*step);
            nn +:= 1;
        end while;
        Append(~units, u);
    end for;
    produ := &*[ SS | units[i]^(r[i]) : i in [1..nd] | r[i] ne 0 ];
    c0 := Coefficient(produ, -L);
    if c0 eq 0 then return CC!0, L, W; end if;
    fac0, z0 := slashdata(w, tau0);
    fac1, z1 := slashdata(w, tau1);
    num0 := fac0^wt * &*[ CC | DedekindEta(d*z0)^(r[i]) : i->d in ds | r[i] ne 0 ];
    num1 := fac1^wt * &*[ CC | DedekindEta(d*z1)^(r[i]) : i->d in ds | r[i] ne 0 ];
    sfun := func< tau | ee(tau*L/(24*W)) *
        &*[ CC | ( DedekindEta((tri[i][1]*tau + tri[i][2])/tri[i][3]) *
                   ee(-(tri[i][1]*tau + tri[i][2])/(24*tri[i][3])) )^(r[i])
            : i in [1..nd] | r[i] ne 0 ] >;
    k0 := num0 / sfun(tau0);
    k1 := num1 / sfun(tau1);
    // RELATIVE consistency check (see eis32k): at large M the eta evaluations
    // lose enough digits that an absolute 1e-30 bound fires spuriously.
    kscale := Abs(k0); if kscale lt 1 then kscale := Parent(kscale)!1; end if;
    error if Abs(k0 - k1) gt 10^(-25)*kscale, "kappa not constant at wt", wt, Abs(k0-k1);
    return k0 * c0, L, W;
end function;

// ---- E constant-term vectors --------------------------------------------
Emat := [ ];  // rows: E's, cols: cosets
for Ei->rE in Epool do
    row := [ CC | 0 : wi in [1..nw] ];
    for wi->w in words do
        v, L, W := a0at(w, rE, 3);
        row[wi] := v;
    end for;
    Append(~Emat, row);
    if Ei mod 10 eq 0 then printf "  E %o/%o done\n", Ei, #Epool; end if;
end for;
printf "E constant table built\n";
// full constants dump for the offline rationalization of the cusp subspace:
// EMAT <E index> <coset index> <class> re im   (50 digits)
R50 := RealField(50);
for Ei in [1..#Epool] do
    for wi->w in words do
        v := Emat[Ei][wi];
        if Abs(v) lt 10^(-40) then continue; end if;
        g := VVWordMatrix(w);
        printf "EMAT %o %o %o %o %o\n", Ei, wi, GCD(g[2][1] mod M, M), R50!Re(v), R50!Im(v);
    end for;
end for;
// rho too, same grid
for wi->w in words do
    if Abs(rhov[wi]) lt 10^(-40) then continue; end if;
    g := VVWordMatrix(w);
    printf "RHOV %o %o %o %o\n", wi, GCD(g[2][1] mod M, M), R50!Re(rhov[wi]), R50!Im(rhov[wi]);
end for;

// ---- least-squares solve: sum_E beta_E a0(E|w) = rho_w -------------------
nE := #Epool;
// no pool supplied: this is a rho-only run (regression / RHOV dump)
if nE eq 0 then printf "no E pool supplied -- rho only\nDONE\n"; quit; end if;
RR := RealField(PREC);
ncol := 2*nE;
col := function(a)  // real 2nw-vector of column a (a <= ne: E_a; else i*E_{a-ne})
    if a le nE then
        return [ RR | Re(Emat[a][wi]) : wi in [1..nw] ] cat [ RR | Im(Emat[a][wi]) : wi in [1..nw] ];
    else
        return [ RR | -Im(Emat[a-nE][wi]) : wi in [1..nw] ] cat [ RR | Re(Emat[a-nE][wi]) : wi in [1..nw] ];
    end if;
end function;
cols := [ col(a) : a in [1..ncol] ];
rvec := [ RR | Re(rhov[wi]) : wi in [1..nw] ] cat [ RR | Im(rhov[wi]) : wi in [1..nw] ];
G2 := ZeroMatrix(RR, ncol, ncol);
for a in [1..ncol] do
    for b in [a..ncol] do
        s := &+[ cols[a][i]*cols[b][i] : i in [1..2*nw] ];
        G2[a][b] := s; G2[b][a] := s;
    end for;
end for;
tr := &+[ G2[a][a] : a in [1..ncol] ];
for a in [1..ncol] do G2[a][a] +:= tr*10^(-45); end for;
cvec := Matrix(RR, 1, ncol, [ &+[ cols[a][i]*rvec[i] : i in [1..2*nw] ] : a in [1..ncol] ]);
ok, X := IsConsistent(G2, cvec);
if ok then
    beta := [ CC!X[1][a] + ii*CC!X[1][nE+a] : a in [1..nE] ];
    resid := Sqrt(&+[ Abs(&+[ beta[a]*Emat[a][wi] : a in [1..nE] ] - rhov[wi])^2 : wi in [1..nw] ]);
    rhonorm := Sqrt(&+[ Abs(rhov[wi])^2 : wi in [1..nw] ]);
    printf "SOLVE resid = %o  (|rho| = %o)\n", RealField(10)!resid, RealField(10)!rhonorm;
    if resid gt 10^(-20)*rhonorm then
        clsres := AssociativeArray(); clsrho := AssociativeArray();
        for wi->w in words do
            g := VVWordMatrix(w);
            cls := GCD(g[2][1] mod M, M);
            dv := Abs(&+[ beta[a]*Emat[a][wi] : a in [1..nE] ] - rhov[wi])^2;
            rv2 := Abs(rhov[wi])^2;
            if IsDefined(clsres, cls) then clsres[cls] +:= dv; clsrho[cls] +:= rv2;
            else clsres[cls] := dv; clsrho[cls] := rv2; end if;
        end for;
        for cls in Sort([ c : c in Keys(clsres) ]) do
            printf "RESCLASS g=%o  |resid|^2 = %o   |rho|^2 = %o\n", cls,
                RealField(10)!clsres[cls], RealField(10)!clsrho[cls];
        end for;
    end if;
    if resid lt 10^(-20)*rhonorm then
        printf "RHO IS IN THE SPAN: E* found\n";
        for a in [1..nE] do
            if Abs(beta[a]) gt 10^(-25) then
                printf "BETA %o %o %o  r=%o\n", a, RealField(30)!Re(beta[a]),
                    RealField(30)!Im(beta[a]), Epool[a];
            end if;
        end for;
        idwi := 0;
        for wi->w in words do
            g := VVWordMatrix(w);
            if g[2][1] mod M eq 0 then idwi := wi; break; end if;
        end for;
        printf "oo-coset index %o\n", idwi;
        QQ := LaurentSeriesRing(CC); q := QQ.1;  // grid q = e(tau/24)
        DEPTH := 24*40;
        eser := AssociativeArray();
        for i->d in ds do
            u := QQ!1 + O(q^DEPTH);
            nn := 1;
            while 24*d*nn lt DEPTH do u *:= 1 - q^(24*d*nn); nn +:= 1; end while;
            eser[d] := q^(d) * u;  // eta(d tau) up to q^{1/24}-grid power d
        end for;
        Estar := QQ!0 + O(q^DEPTH);
        for a in [1..nE] do
            if Abs(beta[a]) lt 10^(-25) then continue; end if;
            rE := Epool[a];
            mono := QQ!1 + O(q^DEPTH);
            for i->d in ds do
                if rE[i] ne 0 then mono *:= eser[d]^rE[i]; end if;
            end for;
            Estar +:= beta[a]*mono;
        end for;
        printf "ESTAR oo-expansion (exponent grid q^(n/24)):\n";
        for nn in [0..DEPTH-1] do
            c := Coefficient(Estar, nn);
            if Abs(c) gt 10^(-25) then
                printf "ECOEF %o/24 %o %o\n", nn, RealField(30)!Re(c), RealField(30)!Im(c);
            end if;
        end for;
    end if;
else
    printf "normal equations singular/inconsistent -- rank deficient pool\n";
end if;
printf "DONE\n";
quit;
