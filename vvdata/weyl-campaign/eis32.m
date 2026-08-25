// eis32: the weight-3/2 side of the residue-pairing route (classfunc-derivation.md).
// Build weight-3/2 eta monomials E with the CONJUGATE panel character as triples
// of panel monomials (Sum r = 3), filter the ones holomorphic at every cusp
// (Ligozat: sum_d r_d gcd(c,d)^2/d >= 0 for all c | M), compute their constant
// terms a0(E|w) at every coset with the cusp6 machinery (weight-general kappa
// pinning), and SOLVE: is rho_w = VVRhoInvE0FFT(w)[eta*] in the span of the
// E-constant vectors?  If yes, print the combination, its oo q-expansion and
// its 0-coset expansion -- W(m) = -a_{E*}(m)/lambda is then the m=0 functional,
// to be checked against the ledger offline.
//
//   magma -b DD:=15 NN:=2 KEYS:=-2,-1,9,10,11,12,13,14,15 eis32.m
AttachSpec("ShimuraQuotients.spec");

D := 15; N := 2; keystr := "-2,-1,9,10,11,12,13,14,15";
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;
if assigned KEYS then keystr := KEYS; end if;
keys := [StringToInteger(s) : s in Split(keystr, ",")];
MAXE := 80; if assigned ME then MAXE := StringToInteger(ME); end if;

PREC := 80;
CC := ComplexField(PREC); ii := CC.1; pi := Pi(CC);
ee := func< z | Exp(2*pi*ii*z) >;

Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
R := Parent(fs[keys[1]]); ds := R`ds;
M := IsOdd(D*N) select 4*D*N else 2*D*N;
printf "BASE %o %o  M = %o  ds = %o\n", D, N, M, ds;

Ld := ShimuraCurveLattice(D, N);
fftdata := VVWeilFFT(Ld, CC : Dual := true);
elts := fftdata[7];
Qr := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
qofidx := function(i)
    v := ChangeRing(elts[i]@@Ld`to_disc, Rationals());
    r := (v*Qr, v)/(2*dn^2); return r - Floor(r);
end function;
isoidx := [ i : i in [1..#elts] | qofidx(i) eq 0 ];
est := isoidx[2];

reps := VVCosetReps(M);
words := [ VVSTWord(g) : g in reps ];
printf "#cosets = %o\n", #reps;

monos := {@ @};
for k in keys do for r in Exponents(fs[k]) do Include(~monos, r); end for; end for;
nm := #monos; nd := #ds;
printf "#monomials = %o\n", nm;

// ---- enumerate holomorphic triples --------------------------------------
cs := Divisors(M);
ordvec := function(r)  // 24*ord at cusp c, up to positive width factors
    return [ &+[ Rationals() | r[i]*GCD(c, ds[i])^2/ds[i] : i in [1..nd] ] : c in cs ];
end function;
mvecs := [ [ Integers() | r[i] : i in [1..nd] ] : r in monos ];
Epool := {@ @};
// The triple enumeration below is CUBIC in #monomials (nm^3/6 ordvec calls):
// fine at nm ~ 200-500, lethal at nm ~ 2000+ (77_2: 1.5e9 iterations, 10+ h).
// Set NOTRIPLES:=1 to skip it and rely entirely on the external EF pool.
skiptriples := assigned NOTRIPLES;
if not skiptriples then
for i in [1..nm] do
    for j in [i..nm] do
        rij := [ mvecs[i][t] + mvecs[j][t] : t in [1..nd] ];
        for k in [j..nm] do
            rE := [ rij[t] + mvecs[k][t] : t in [1..nd] ];
            if forall{ o : o in ordvec(rE) | o ge 0 } then
                Include(~Epool, rE);
                if #Epool ge MAXE then break i; end if;
            end if;
        end for;
    end for;
end for;
end if;
printf "#holomorphic triple E's = %o (cap %o)\n", #Epool, MAXE;
// optional external pool (python-enumerated holomorphic quotients): EF := path
// to a file containing a Magma sequence literal [ [r_1,...,r_nd], ... ]
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
error if #Epool eq 0, "no holomorphic E's -- eta span insufficient, need honest Eisenstein";

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

// a0 of the eta monomial with exponent vector r (weight (sum r)/2) at coset w;
// if depth0 > 0 also print the first expansion coefficients (for E* readout).
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

// ---- E constant-term vectors + rho vector -------------------------------
nw := #words;
rhov := [ CC | 0 : wi in [1..nw] ];
for wi->w in words do
    rv := VVRhoInvE0FFT(fftdata, w);
    rhov[wi] := rv[est];
end for;

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
