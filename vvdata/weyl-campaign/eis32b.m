// eis32b: END-TO-END acceptance of the weight-3/2 residue pairing on a base.
// Re-solves E* (as eis32, from the EF pool), then computes, for every panel
// form f, A_f = sum_w a0((f * E*)|w) over ALL cosets -- the residue theorem
// demands A_f = 0 with every convention (multiplier conjugacy, weight-2 slash,
// fractional grids) live.  Also prints the block decomposition
//   INTERM_f = sum over cosets where every panel monomial is holomorphic,
//   CONV_f   = the convolution blocks (cosets with monomial poles: oo/0 class)
// CORRECTION (2026-08-28): the claim that used to stand here -- "so INTERM_f =
// c_eta*(0) and CONV_f = -c_eta*(0) exhibit the functional" -- is WRONG, and
// was wrong in every run ever made with this script.  The INTERM/CONV split is
// by PANEL-MONOMIAL pole support, a property of the panel and not of f: a coset
// can carry a monomial pole while f itself is holomorphic there.  Every log
// disagrees with the old claim (15_2 FORM -2: INTERM ~ 0 against c_eta*(0) = 4;
// FORM 9: INTERM = -8 against 0).  The identity that DOES hold is the
// n = 0 / n > 0 split of the residue sum,
//   c_eta*(0) = sum_w a_0(f|w) a_0(E|w) = - sum_w sum_{n>0} a_{-n}(f|w) a_n(E|w),
// proved as "the residue pairing" in the preprint.  A_f = 0 is the residue
// theorem itself and is untouched by this correction.
// Finally dumps E*'s expansion at one canonical coset per class (E0COEF lines)
// -- the 0-class lines are the ledger's 0-side weights.
//
//   magma -b DD:=15 NN:=2 EF:=vvdata/weyl-campaign/epool_15_2.txt eis32b.m
AttachSpec("ShimuraQuotients.spec");

D := 15; N := 2; keystr := "-2,-1,9,10,11,12,13,14,15";
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;
if assigned KEYS then keystr := KEYS; end if;
keys := [StringToInteger(s) : s in Split(keystr, ",")];

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
nd := #ds;
printf "BASE %o %o  M = %o\n", D, N, M;

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
nw := #words;

monos := {@ @};
for k in keys do for r in Exponents(fs[k]) do Include(~monos, r); end for; end for;
nm := #monos;
mvecs := [ [ Integers() | r[i] : i in [1..nd] ] : r in monos ];

error if not assigned EF, "EF (E-pool file) required";
Epool := {@ @};
for r in eval Read(EF) do Include(~Epool, [ Integers() | x : x in r ]); end for;
nE := #Epool;
printf "#panel monomials = %o, #E pool = %o\n", nm, nE;

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
    return a, b, e;
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

// per-coset machinery, split so pass 1 never builds deep series it won't use.
// cosetleads: tri + leads only (cheap).  cosetfull: units once, then produ/kappa
// per vector only where need[vi] is set; needk[vi] requests kappa alone
// (produ replaced by 1 -- valid ONLY for reading the t^0 coefficient at L = 0).
cosetleads := function(w, vecs)
    g := VVWordMatrix(w);
    tri := [ ];
    for d in ds do
        a, b, e := triang(g, d);
        Append(~tri, <a, b, e>);
    end for;
    W := LCM([ tri[i][3] : i in [1..nd] ]);
    leads := [ &+[ Integers() | r[i]*tri[i][1]*(W div tri[i][3]) : i in [1..nd] ] : r in vecs ];
    return leads, tri, W;
end function;

cosetfull := function(w, vecs, tri, W, leads, need, needk)
    depth := Maximum([ 1 ] cat [ -leads[vi] + 1 : vi in [1..#vecs] | need[vi] ]);
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
    fac0, z0 := slashdata(w, tau0);
    fac1, z1 := slashdata(w, tau1);
    sers := [ SS | ]; kaps := [ CC | ];
    for vi->r in vecs do
        if not (need[vi] or needk[vi]) then
            Append(~sers, SS!1); Append(~kaps, CC!0); continue;
        end if;
        L := leads[vi];
        produ := need[vi]
            select &*[ SS | units[i]^(r[i]) : i in [1..nd] | r[i] ne 0 ]
            else SS!1;
        wt := &+[ Integers() | x : x in r ];
        num0 := fac0^wt * &*[ CC | DedekindEta(d*z0)^(r[i]) : i->d in ds | r[i] ne 0 ];
        num1 := fac1^wt * &*[ CC | DedekindEta(d*z1)^(r[i]) : i->d in ds | r[i] ne 0 ];
        sfun := func< tau | ee(tau*L/(24*W)) *
            &*[ CC | ( DedekindEta((tri[i][1]*tau + tri[i][2])/tri[i][3]) *
                       ee(-(tri[i][1]*tau + tri[i][2])/(24*tri[i][3])) )^(r[i])
                : i in [1..nd] | r[i] ne 0 ] >;
        k0 := num0 / sfun(tau0);
        k1 := num1 / sfun(tau1);
        error if Abs(k0 - k1) gt 10^(-30), "kappa not constant";
        Append(~sers, produ); Append(~kaps, k0);
    end for;
    return sers, kaps;
end function;

allvecs := mvecs cat [ r : r in Epool ];   // panel monomials then E's
Af    := [ CC | 0 : k in keys ];   // total residue sum per form
Aint  := [ CC | 0 : k in keys ];   // holomorphic-coset (factorized) block
Aconv := [ CC | 0 : k in keys ];   // convolution block (cosets w/ poles)
ceta  := [ CC | 0 : k in keys ];   // assembly check sum rho a0(f|w)

// solve beta first: need a0(E|w) for all cosets -- collect while streaming
Emat := [ [ CC | 0 : wi in [1..nw] ] : a in [1..nE] ];
f0tab := [ [ CC | 0 : wi in [1..nw] ] : ri in [1..nm] ];  // a0(m_r|w)
convterm := [ [ CC | 0 : wi in [1..nw] ] : ri in [1..nm] ]; // a0(m_r E*|w) filled pass 2
Ldata := [ ]; // store leads per coset for pass 2 reuse decision
printf "pass 1: constants\n";
for wi->w in words do
    leads, tri, W := cosetleads(w, allvecs);
    need  := [ vi le nm and leads[vi] lt 0 : vi in [1..#allvecs] ];
    needk := [ leads[vi] eq 0 : vi in [1..#allvecs] ];
    sers, kaps := cosetfull(w, allvecs, tri, W, leads, need, needk);
    for ri in [1..nm] do
        L := leads[ri];
        if L gt 0 then f0tab[ri][wi] := 0;
        elif L eq 0 then f0tab[ri][wi] := kaps[ri];
        else f0tab[ri][wi] := kaps[ri]*Coefficient(sers[ri], -L); end if;
    end for;
    for a in [1..nE] do
        L := leads[nm+a];
        error if L lt 0, "E pool member not holomorphic at coset", wi;
        Emat[a][wi] := L eq 0 select kaps[nm+a] else CC!0;
    end for;
    Append(~Ldata, leads);
    if wi mod 24 eq 0 then printf "  coset %o/%o\n", wi, nw; end if;
end for;

rhov := [ CC | 0 : wi in [1..nw] ];
for wi->w in words do
    rv := VVRhoInvE0FFT(fftdata, w);
    rhov[wi] := rv[est];
end for;

// realified least squares for beta (as eis32)
RR := RealField(PREC);
ncol := 2*nE;
col := function(a)
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
error if not ok, "solve failed";
beta := [ CC!X[1][a] + ii*CC!X[1][nE+a] : a in [1..nE] ];
resid := Sqrt(&+[ Abs(&+[ beta[a]*Emat[a][wi] : a in [1..nE] ] - rhov[wi])^2 : wi in [1..nw] ]);
printf "SOLVE resid = %o\n", RealField(10)!resid;

// pass 2: per coset, E* series and convolutions where panel monomials have poles
printf "pass 2: convolutions\n";
e0dumped := {};
for wi->w in words do
    leads := Ldata[wi];
    anypole := exists{ri : ri in [1..nm] | leads[ri] lt 0};
    g := VVWordMatrix(w);
    cls := GCD(g[2][1] mod M, M);
    if not anypole then
        // factorized: a0(m_r E*) = a0(m_r) * a0(E*)
        estar0 := &+[ CC | beta[a]*Emat[a][wi] : a in [1..nE] ];
        for ri in [1..nm] do
            convterm[ri][wi] := f0tab[ri][wi] * estar0;
        end for;
        continue;
    end if;
    depth := Maximum([ -leads[ri] : ri in [1..nm] ]) + 1;
    leads2, tri, W := cosetleads(w, allvecs);
    need  := [ (vi le nm and leads2[vi] lt 0) or (vi gt nm and leads2[vi] lt depth) : vi in [1..#allvecs] ];
    // panel monomials with L = 0 still contribute kappa * Eser(0): kappa needed
    needk := [ vi le nm and leads2[vi] eq 0 : vi in [1..#allvecs] ];
    sers, kaps := cosetfull(w, allvecs, tri, W, leads2, need, needk);
    Eser := SS!0 + O(t^depth);
    for a in [1..nE] do
        LE := leads2[nm+a];
        if LE ge depth then continue; end if;
        Eser +:= beta[a]*kaps[nm+a]*(t^LE*sers[nm+a] + O(t^depth));
    end for;
    // dump E*'s expansion once per class (ledger 0-side weights)
    if cls notin e0dumped then
        Include(~e0dumped, cls);
        R30 := RealField(30);
        for j in [0..Minimum(depth-1, 24*W)] do
            cj := Coefficient(Eser, j);
            if Abs(cj) gt 10^(-25) then
                printf "E0COEF g=%o W=%o %o/%o %o %o\n", cls, W, j, 24*W, R30!Re(cj), R30!Im(cj);
            end if;
        end for;
    end if;
    for ri in [1..nm] do
        L := leads2[ri];
        if L gt 0 then convterm[ri][wi] := 0; continue; end if;
        // a0(m_r E* | w) = kap_r * sum_j c_j(produ_r) * e_{-L-j}(Eser)
        s := CC!0;
        for j in [0..-L] do
            cj := Coefficient(sers[ri], j);
            if cj ne 0 then s +:= cj*Coefficient(Eser, -L-j); end if;
        end for;
        convterm[ri][wi] := kaps[ri]*s;
    end for;
    printf "  conv coset %o/%o class %o done\n", wi, nw, cls;
end for;

// assemble per form
for ki->k in keys do
    f := fs[k];
    for wi in [1..nw] do
        leads := Ldata[wi];
        tot := CC!0; itot := CC!0;
        for r in Exponents(f) do
            ri := Index(monos, r);
            tot +:= f`coeffs[r]*convterm[ri][wi];
        end for;
        Af[ki] +:= tot;
        anypole := exists{ri : ri in [1..nm] | leads[ri] lt 0};
        if anypole then Aconv[ki] +:= tot; else Aint[ki] +:= tot; end if;
        for r in Exponents(f) do
            ri := Index(monos, r);
            ceta[ki] +:= f`coeffs[r]*rhov[wi]*f0tab[ri][wi];
        end for;
    end for;
    printf "FORM %o  A_f = (%o, %o)  INTERM = (%o, %o)  CONV = (%o, %o)  c_eta*(0) = (%o, %o)\n",
        k, RealField(15)!Re(Af[ki]), RealField(15)!Im(Af[ki]),
        RealField(15)!Re(Aint[ki]), RealField(15)!Im(Aint[ki]),
        RealField(15)!Re(Aconv[ki]), RealField(15)!Im(Aconv[ki]),
        RealField(15)!Re(ceta[ki]), RealField(15)!Im(ceta[ki]);
end for;
printf "DONE\n";
quit;
