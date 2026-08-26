// eis32k: eis32fast with the rhov loop replaced by the KUBOTA CLOSED FORM
// (rhokubota8/9, banked b130990) -- no FFT, no (DN)^3 wall:
//   rv[c y + x_c] = u(w) xi(a,c) sqrt(|D[c]|/|D|) e(ac Q(y) + B(x_c, y)),
//   gi = (word matrix)^{-1} = (a,b;c,d), u = Kubota cocycle sign,
//   x_c per Stromberg Def 2.15 (Jordan-basis sum of the 2-part at v2(c)=1,
//   0 otherwise for an elementary (Z/2)^r 2-part), xi(a,c) per Def 6.1.
// The tracked coefficient is rv[est]: need y with c y = elts[est] - x_c;
// no solution => 0.  c = 0 words (only T^0 among Gamma_0 reps) => 1.
// Requires an ELEMENTARY (Z/2)^3 2-part of the discriminant group (errors
// out otherwise -- extend Def 2.15 handling before odd-D/4-component bases).
//   magma -b DD:=330 NN:=1 EF:=<pool> eis32k.m
AttachSpec("ShimuraQuotients.spec");

D := 15; N := 2;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;

PREC := 80;
CC := ComplexField(PREC); ii := CC.1; pi := Pi(CC);
ee := func< z | Exp(2*pi*ii*z) >;

M := IsOdd(D*N) select 4*D*N else 2*D*N;
ds := Divisors(M); nd := #ds;
printf "BASE %o %o  M = %o  ds = %o\n", D, N, M, ds;
Ld := ShimuraCurveLattice(D, N);
fftdata := VVWeilFFT(Ld, CC : Dual := true);
elts := fftdata[7]; i0 := fftdata[8];
n := #elts;
printf "|D| = %o  i0 = %o\n", n, i0;

Qr := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
lift := [ ChangeRing(elts[i]@@Ld`to_disc, Rationals()) : i in [1..n] ];
frac := func< r | r - Floor(r) >;
Qbar := [ frac(-(lift[i]*Qr, lift[i])/(2*dn^2)) : i in [1..n] ];
BB := function(i, j)
    return frac(-(lift[i]*Qr, lift[j])/dn^2);
end function;
idx := AssociativeArray();
for i in [1..n] do idx[elts[i]] := i; end for;

isoidx := [ i : i in [1..n] | Qbar[i] eq 0 ];
// #iso = 2N-1: at N=1 the zero coset is the ONLY isotropic element, so track
// it; for N>1 keep the first nonzero isotropic coset (unchanged behavior).
est := #isoidx ge 2 select isoidx[2] else isoidx[1];
printf "#isotropic = %o, tracking coset index %o\n", #isoidx, est;

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

// ---- rho vector via the closed form ------------------------------------
nw := #words;
rhov := [ CC | 0 : wi in [1..nw] ];
// per-c-residue: a preimage y0 with c*y0 = elts[est] - x_c (or 0 = no
// support at est), and |D[c]|
YTab := AssociativeArray();  // <cm, xci> -> <y0 index or 0, dcn>
for wi->w in words do
    if #w eq 0 or VVWordMatrix(w)[2][1] eq 0 then
        // T^k coset: rho^{-1} diagonal, e(k Q(est)) = 1 (est isotropic)
        rhov[wi] := CC!1;
        continue;
    end if;
    gmat := VVWordMatrix(w);
    gi := gmat^(-1);
    a := gi[1][1]; b := gi[1][2]; c := gi[2][1]; d := gi[2][2];
    if a eq 0 then
        // the S-coset: the closed form carries ONE base-dependent global
        // sign here (-1 at 15_2/10_7, +1 at 210_1) -- measure it by FFT
        // (a single word, affordable) instead of guessing the bit
        rv := VVRhoInvE0FFT(fftdata, w);
        rhov[wi] := rv[est];
        printf "AZERO wi=%o rhov(FFT) = %o + %o i\n", wi,
            RealField(30)!Re(rv[est]), RealField(30)!Im(rv[est]);
        continue;
    end if;
    k2 := Valuation(c, 2);
    xci := (k2 eq 1) select xc2i else i0;
    cm := ((c mod levD) + levD) mod levD;
    key := <cm, xci>;
    if not IsDefined(YTab, key) then
        tgt := elts[est] - (xci eq i0 select Parent(elts[1])!0 else elts[xci]);
        y0 := 0;
        dcn := 0;
        for j in [1..n] do
            if IsZero(cm*elts[j]) then dcn +:= 1; end if;
            if y0 eq 0 and cm*elts[j] eq tgt then y0 := j; end if;
        end for;
        YTab[key] := <y0, dcn>;
    end if;
    y0 := YTab[key][1]; dcn := YTab[key][2];
    if y0 eq 0 then continue; end if;   // rhov[wi] stays 0
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
    rhov[wi] := u * xidef * Sqrt(CC!dcn/CC!n) *
            ee(frac(a*c*Qbar[y0] + (xci eq i0 select 0 else BB(xci, y0))));
    if wi mod 128 eq 0 then printf "rho %o/%o (closed form)\n", wi, nw; end if;
end for;
printf "rhov done (closed form, %o YTab keys)\n", #Keys(YTab);

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
        n := 1;
        while n*step lt depth do
            u *:= 1 - ee(CC!(n*b/e))*t^(n*step);
            n +:= 1;
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
    error if Abs(k0 - k1) gt 10^(-30), "kappa not constant at wt", wt, Abs(k0-k1);
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
// realified least squares (beta = u + i v): columns [Re E; Im E] and
// [-Im E; Re E]; real symmetric Gram + small ridge (pool is expected to be
// rank-deficient); the residual is the arbiter, so the ridge bias is safe.
nE := #Epool;
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
G := ZeroMatrix(RR, ncol, ncol);
for a in [1..ncol] do
    for b in [a..ncol] do
        s := &+[ cols[a][i]*cols[b][i] : i in [1..2*nw] ];
        G[a][b] := s; G[b][a] := s;
    end for;
end for;
tr := &+[ G[a][a] : a in [1..ncol] ];
for a in [1..ncol] do G[a][a] +:= tr*10^(-45); end for;
cvec := Matrix(RR, 1, ncol, [ &+[ cols[a][i]*rvec[i] : i in [1..2*nw] ] : a in [1..ncol] ]);
ok, X := IsConsistent(G, cvec);
if ok then
    beta := [ CC!X[1][a] + ii*CC!X[1][nE+a] : a in [1..nE] ];
    resid := Sqrt(&+[ Abs(&+[ beta[a]*Emat[a][wi] : a in [1..nE] ] - rhov[wi])^2 : wi in [1..nw] ]);
    rhonorm := Sqrt(&+[ Abs(rhov[wi])^2 : wi in [1..nw] ]);
    printf "SOLVE resid = %o  (|rho| = %o)\n", RealField(10)!resid, RealField(10)!rhonorm;
    if resid gt 10^(-20)*rhonorm then
        // where does the miss live?  residual and rho norms by cusp class
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
        // E* oo-expansion: identity coset is the one with word []; find it
        idwi := 0;
        for wi->w in words do
            g := VVWordMatrix(w);
            if g[2][1] mod M eq 0 then idwi := wi; break; end if;
        end for;
        printf "oo-coset index %o\n", idwi;
        // oo q-expansion of E* = sum beta_E eta-monomial(rE): exact eta series
        QQ := LaurentSeriesRing(CC); q := QQ.1;  // grid q = e(tau/24)
        DEPTH := 24*40;
        eser := AssociativeArray();
        for i->d in ds do
            u := QQ!1 + O(q^DEPTH);
            n := 1;
            while 24*d*n lt DEPTH do u *:= 1 - q^(24*d*n); n +:= 1; end while;
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
        for n in [0..DEPTH-1] do
            c := Coefficient(Estar, n);
            if Abs(c) gt 10^(-25) then
                printf "ECOEF %o/24 %o %o\n", n, RealField(30)!Re(c), RealField(30)!Im(c);
            end if;
        end for;
    end if;
else
    printf "normal equations singular/inconsistent -- rank deficient pool\n";
end if;
printf "DONE\n";
quit;
