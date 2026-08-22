// Cusp-expansion machinery for the m = 0 multiplier: c_{eta*}(0) as a FINITE sum
// over Gamma_0(M)\SL_2(Z) cosets, replacing the numeric oracle's Fourier sampling.
//
// For each coset word w (matrix g) and eta monomial m_r = prod_d eta(d tau)^{r_d}:
//   [d 0; 0 1] g = g_d * [a_d b_d; 0 e_d]   (g_d in SL_2(Z), a_d e_d = d)
//   (m_r | w)(tau) = kappa * Q^{alpha} prod_d [ prod_n (1 - e(n b_d/e_d) Q^{n a_d/e_d}) ]^{r_d}
// with Q = e(tau), alpha = sum_d r_d a_d/(24 e_d), and kappa a CONSTANT pinned
// numerically at one point (checked at a second) -- the eta-multiplier/metaplectic
// bookkeeping never enters.  The q^0 coefficient of f|w is then exact-in-structure,
// and c_{eta*}(0) = sum_w (rho(w^-1)e_0)_{eta*} a_0(f|w).   mult = c_{eta*}(0)/2.
//
//   magma -b DD:=15 NN:=2 KEYS:=-1,13 cusp1.m
AttachSpec("ShimuraQuotients.spec");

D := 15; N := 2; keystr := "-1,13";
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;
if assigned KEYS then keystr := KEYS; end if;
keys := [StringToInteger(s) : s in Split(keystr, ",")];

PREC := 80;
CC := ComplexField(PREC); ii := CC.1; pi := Pi(CC);
ee := func< z | Exp(2*pi*ii*z) >;   // e(z) for complex z

// ---- the forms ----
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
elts := fftdata[7]; i0 := fftdata[8];
// isotropic cosets of A
Qr := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
qofidx := function(i)
    v := ChangeRing(elts[i]@@Ld`to_disc, Rationals());
    r := (v*Qr, v)/(2*dn^2); return r - Floor(r);
end function;
isoidx := [ i : i in [1..#elts] | qofidx(i) eq 0 ];
printf "#isotropic = %o (2N-1 = %o)\n", #isoidx, 2*N-1;

reps := VVCosetReps(M);
words := [ VVSTWord(g) : g in reps ];
printf "#cosets = %o\n", #reps;

// ---- per-coset triangularization data ----
// for g and each d: g_d, (a_d, b_d, e_d) with [d 0;0 1] g = g_d [a b;0 e]
triang := function(g, d)
    g2 := Matrix(Integers(), 2, 2, [d*g[1][1], d*g[1][2], g[2][1], g[2][2]]);
    c1 := g2[1][1]; c2 := g2[2][1];
    h := GCD(c1, c2);
    if h eq 0 then error "degenerate"; end if;
    p1 := c1 div h; p2 := c2 div h;
    _, y, x := XGCD(p1, -p2);          // p1*y - p2*x ... solve p1*y2 - p2*x2 = 1
    // find (x2,y2) with p1*y2 - x2*p2 = 1: XGCD(p1, p2) = <1, u, v>: u p1 + v p2 = 1
    gg, u, v := XGCD(p1, p2);
    error if gg ne 1, "not primitive";
    gd := Matrix(Integers(), 2, 2, [p1, -v, p2, u]);   // det = p1 u + v p2 = 1
    sd := gd^(-1) * g2;                                 // upper triangular, [a b; 0 e]
    error if sd[2][1] ne 0 or sd[1][1]*sd[2][2] ne d, "triangularization failed";
    a := sd[1][1]; b := sd[1][2]; e := sd[2][2];
    if a lt 0 then a := -a; b := -b; e := -e; end if;   // normalize a > 0
    if e lt 0 then                                       // want e > 0: flip gd sign
        gd := -gd; sd := -sd; a := sd[1][1]; b := sd[1][2]; e := sd[2][2];
    end if;
    return a, b, e, gd;
end function;

// numeric slash factor of a word at tau (phi_w(tau)^{-1}), and the moved point
slashdata := function(word, tau)
    z := tau; factor := CC!1;
    for i := #word to 1 by -1 do
        if word[i][1] eq "S" then factor /:= Sqrt(z); z := -1/z;
        else z := z + word[i][2]; end if;
    end for;
    return factor, z;   // (f|w)(tau) = factor * f(z), z = g(tau)
end function;

tau0 := CC!0.31 + CC!1.31*ii;
tau1 := CC!(-0.57) + CC!1.73*ii;

// union of exponent vectors over requested forms
monos := {@ @};
for k in keys do for r in Exponents(fs[k]) do Include(~monos, r); end for; end for;
printf "#monomials = %o\n", #monos;

// ---- per coset: compute a0 data for every monomial ----
// store A0[w][r-index] = kappa * (coefficient of Q^0), as a complex number
SS := PowerSeriesRing(CC); t := SS.1;
a0tab := [ [ CC!0 : r in monos ] : w in words ];
for wi->w in words do
    g := VVWordMatrix(w);
    tri := [ ]; // per d
    for d in ds do
        a, b, e, gd := triang(g, d);
        Append(~tri, <a, b, e, gd>);
    end for;
    W := LCM([ tri[i][3] : i in [1..#ds] ]);
    // leading t-exponents per monomial (t = Q^{1/24W})
    leads := [ &+[ Integers() | r[i]*tri[i][1]*(W div tri[i][3]) : i in [1..#ds] ]
               : r in monos ];
    depth := Maximum([ 0 ] cat [ -L : L in leads ]) + 1;
    // unit series per d to precision depth
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
    for ri->r in monos do
        L := leads[ri];
        // A0 = coefficient of t^{-L}... constant term of t^L prod u^r: coeff of t^0
        // i.e. coefficient of t^{-L} in prod_d u_d^{r_d}; zero when L > 0 or off grid
        if L gt 0 then continue; end if;
        produ := &*[ SS | units[i]^(r[i]) : i in [1..#ds] | r[i] ne 0 ];
        c0 := Coefficient(produ, -L);
        if c0 eq 0 and -L gt 0 then a0tab[wi][ri] := CC!0; continue; end if;
        // kappa: numeric ratio at tau0, verified at tau1
        num0 := fac0 * &*[ CC | DedekindEta(d*z0)^(r[i]) : i->d in ds | r[i] ne 0 ];
        num1 := fac1 * &*[ CC | DedekindEta(d*z1)^(r[i]) : i->d in ds | r[i] ne 0 ];
        sfun := func< tau | ee(tau*L/(24*W)) *
            &*[ CC | ( DedekindEta((tri[i][1]*tau + tri[i][2])/tri[i][3]) *
                       ee(-(tri[i][1]*tau + tri[i][2])/(24*tri[i][3])) )^(r[i])
                : i in [1..#ds] | r[i] ne 0 ] >;
        k0 := num0 / sfun(tau0);
        k1 := num1 / sfun(tau1);
        error if Abs(k0 - k1) gt 10^(-PREC+25), "kappa not constant", wi, ri, Abs(k0-k1);
        a0tab[wi][ri] := k0 * c0;
    end for;
end for;
printf "a0 table built\n";

// ---- assemble c_eta(0) per form ----
for k in keys do
    f := fs[k];
    a0w := [ &+[ CC | f`coeffs[r] * a0tab[wi][Index(monos, r)] : r in Exponents(f) ]
             : wi in [1..#words] ];
    // c_eta(0) = sum_w rhovec_w[eta] * a0(f|w)
    cvals := [ CC!0 : i in [1..#isoidx] ];
    for wi->w in words do
        if Abs(a0w[wi]) lt 10^(-PREC+10) then continue; end if;
        rv := VVRhoInvE0FFT(fftdata, w);
        for j->i in isoidx do cvals[j] +:= rv[i] * a0w[wi]; end for;
    end for;
    printf "FORM %o:\n", k;
    for j->i in isoidx do
        tag := i eq i0 select " (trivial)" else "";
        printf "  c_eta(0)[%o]%o = %o + %o i   (half: %o)\n", j, tag,
            RealField(12)!Re(cvals[j]), RealField(12)!Im(cvals[j]),
            RealField(12)!(Re(cvals[j])/2);
    end for;
end for;
quit;
