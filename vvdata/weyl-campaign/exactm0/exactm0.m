// The m = 0 multipliers (1/2) c_eta(0) of Schofer's formula computed with NO transcendental step.
//
// c_eta(0) = sum over the cosets w of Gamma_0(M)\SL_2(Z) of  rho*(w^-1) e_0 [eta] * a_0(f | w),  and the
// per-coset product is constant on the class g = gcd(c, M) (thetag-derivation.md; what
// M0MultipliersBySupport samples numerically).  Here each class is evaluated from ONE representative
// gamma = [a b; c d] with c > 0 and the canonical lift (gamma, sqrt(c tau + d)), with every factor
// an exact algebraic number:
//
//   rho*(gamma^-1) e_0 [eta] = rho_B(gamma)_{0, eta}
//        = e(1/8) c^{-3/2} |L^v/L|^{-1/2} e(d Q(eta)/c) G(c, a; eta),
//        G(c, a; eta) = sum_{nu in L/cL} e((a Q(nu) + (eta, nu))/c),            [theta transformation]
//   f | gamma = e(-1/8) sum_r c_r prod_d [eps(g_d) e(b_d/(24 e_d)) e_d^{-1/2}]^{r_d}
//                                    prod_d [e(a_d tau/(24 e_d)) prod_n (1 - e(n b_d/e_d) q^{n a_d/e_d})]^{r_d}
//        with [d 0; 0 1] gamma = g_d [a_d b_d; 0 e_d], eps the Dedekind-eta multiplier of Apostol,
//        Thm 3.4 (eta(g tau) = eps(g) (-i(c tau + d))^{1/2} eta(tau), c > 0), a 24th root of unity
//        given by a Dedekind sum.  (cusp4.m checked this closed form against the numerically pinned
//        slash constants of M0MultipliersBySupport on all 144 cosets of X_0^15(2).)
//
// Roots of unity live in mu_n, n = 24 M, and are handled as monomials of Q[x]/(x^n - 1); the square
// roots of rationals are the positive real roots, written as Gauss sums in the same ring; the result is
// reduced modulo the n-th cyclotomic polynomial and must be a rational number.  Verify := true
// recomputes every class 1 < g < M from a second representative with a different bottom row and
// demands exact agreement -- an internal check of both closed forms.
//
//   magma -b DD:=15 NN:=2 exactm0.m                 (compare with M0MultipliersBySupport; ~10 s)
//   magma -b DD:=15 NN:=2 RHOCHECK:=1 exactm0.m     (also: the rho formula against the FFT, numerically)
//   magma -b DD:=6 NN:=35 CMP:=0 exactm0.m          (exact only)
//   magma -b DD:=6 NN:=35 CMP:=0 LLLFORMS:=1 exactm0.m   (the four compmult5 probe forms instead of BorcherdsForms)
// Run it from a checkout that has M0MultipliersBySupport (branch composite-level) when CMP is on.
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);

D := 15; N := 2; cmp := true; rhocheck := false; verify := true; lll := false;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;
if assigned CMP then cmp := CMP eq "1"; end if;
if assigned RHOCHECK then rhocheck := RHOCHECK eq "1"; end if;
if assigned VERIFY then verify := VERIFY eq "1"; end if;
if assigned LLLFORMS then lll := LLLFORMS eq "1"; end if;      // (LLL itself is the intrinsic)

// ---------------------------------------------------------------------------------------------
// exact ingredients
// ---------------------------------------------------------------------------------------------

// Dedekind sum s(h, k), exact rational
function dedsum(h, k)
    s := Rationals()!0;
    for i := 1 to k - 1 do
        x := Rationals()!i/k;
        y := Rationals()!(h*i)/k; y := y - Floor(y);
        if y ne 0 then s +:= (x - 1/2)*(y - 1/2); end if;
    end for;
    return s;
end function;

// eps(g) = e(x), x in (1/24)Z, for g in SL_2(Z) with c > 0 (Apostol Thm 3.4); returns x in [0, 1)
function epsexp(g)
    a := g[1][1]; c := g[2][1]; d := g[2][2];
    error if c le 0, "epsexp needs c > 0";
    x := (Rationals()!(a + d)/(12*c) + dedsum(-d, c))/2;
    x := x - Floor(x);
    error if not IsIntegral(24*x), "the eta multiplier is not a 24th root of unity";
    return x;
end function;

// [d 0; 0 1] g = gd * [a b; 0 e] with gd in SL_2(Z), gd[2][1] > 0, a = gcd(d a_g, c) > 0, e = d/a > 0
function triang(g, d)
    c := g[2][1];
    error if c le 0, "triang needs c > 0";
    g2 := Matrix(Integers(), 2, 2, [d*g[1][1], d*g[1][2], c, g[2][2]]);
    h := GCD(g2[1][1], c);
    p1 := g2[1][1] div h; p2 := c div h;
    _, u, v := XGCD(p1, p2);                         // p1 u + p2 v = 1
    gd := Matrix(Integers(), 2, 2, [p1, -v, p2, u]);
    sd := gd^-1 * g2;
    assert Determinant(gd) eq 1 and sd[2][1] eq 0 and sd[1][1] eq h and sd[2][2] eq d div h and gd*sd eq g2;
    return sd[1][1], sd[1][2], sd[2][2], gd;
end function;

// The coefficient of t^(-L) in prod_d P(zeta_{e_d}^{b_d} t^{step_d})^{r_d}, t = q^{1/(24 W)},
// step_d = 24 a_d W/e_d, P(y) = prod_{n >= 1} (1 - y^n) (Euler's pentagonal series), as an element of
// Q(zeta_{W'}) (W' = lcm of the orders of the e(b_d/e_d)), together with that field and W'.
function unit_constant(r, tri, W, L)
    depth := -L + 1;
    act := [i : i in [1..#tri] | r[i] ne 0];
    ords := [ tri[i][3] div GCD(tri[i][2] mod tri[i][3], tri[i][3]) : i in act ];
    Wp := LCM(ords);
    if Wp le 2 then
        F := Rationals(); zW := [F | (-1)^j : j in [0..Wp-1]];
    else
        F := CyclotomicField(Wp); zW := [F.1^j : j in [0..Wp-1]];
    end if;
    S := PowerSeriesRing(F, depth);
    prod := S!1;
    for i in act do
        a, b, e := Explode(tri[i]);
        step := 24*a*(W div e);
        ydepth := (depth - 1) div step;
        Z := PowerSeriesRing(Integers(), ydepth + 1); y := Z.1;
        Pe := Z!1; k := 1;
        while k*(3*k - 1) div 2 le ydepth do
            Pe +:= (-1)^k * y^(k*(3*k - 1) div 2);
            if k*(3*k + 1) div 2 le ydepth then Pe +:= (-1)^k * y^(k*(3*k + 1) div 2); end if;
            k +:= 1;
        end while;
        Pr := r[i] ge 0 select Pe^r[i] else (1/Pe)^(-r[i]);
        bb := b mod e; gg := GCD(bb, e); ord := e div gg; kk := (bb div gg) * (Wp div ord);   // e(b/e) = zeta_W'^kk
        coeffs := [F | 0 : j in [1..depth]];
        for j in [0..ydepth] do
            cj := Coefficient(Pr, j);
            if cj ne 0 then coeffs[j*step + 1] := cj * zW[((kk*j) mod Wp) + 1]; end if;
        end for;
        prod *:= S!coeffs;
    end for;
    return Coefficient(prod, -L), F, Wp;
end function;

// ---------------------------------------------------------------------------------------------
// the exact multiplier
// ---------------------------------------------------------------------------------------------

// returns: per form an associative array p -> (1/2) c_eta(0) for eta of order p (p | N), and the
// per-support-class values (including the mixed classes and the zero coset) for printing
function exact_m0(fs, Ld, D, N : Verify := true, EtaPerClass := 3)
    R := Parent(fs[1]); ds := R`ds; M := R`M;
    assert M eq (IsOdd(D*N) select 4*D*N else 2*D*N);
    n := 24*M;
    P<x> := PolynomialRing(Rationals());
    Phi := CyclotomicPolynomial(n);
    fold := function(h)
        if Degree(h) lt n then return h; end if;
        cs := Coefficients(h); out := [Rationals() | 0 : i in [1..n]];
        for i->ci in cs do out[((i - 1) mod n) + 1] +:= ci; end for;
        return P!out;
    end function;
    zet := func< k | x^(k mod n) >;                                  // e(k/n)
    // the positive real sqrt of a squarefree m >= 1 as a Gauss sum (needs 8 | n and p | n for p | m)
    sqrtpoly := function(m)
        s := P!1;
        if IsEven(m) then s := fold(s*(zet(n div 8) + zet(-(n div 8)))); m := m div 2; end if;
        for p in PrimeDivisors(m) do
            error if n mod p ne 0, "sqrt(p) is not in Q(zeta_n)";
            g := &+[ P | KroneckerSymbol(a, p)*zet(a*(n div p)) : a in [1..p-1] ];
            if p mod 4 eq 3 then g := fold(g*zet(-(n div 4))); end if;
            s := fold(s*g);
        end for;
        return s;
    end function;

    // the discriminant group, lifts v (scaled by dn) with Q(eta) = (v Q v)/(2 dn^2), (eta, nu) = (v Q nu)/dn
    Qm := ChangeRing(Ld`Q, Integers()); dn := Ld`denom;
    elts := [g : g in Ld`disc_grp];
    vs := [];
    for e in elts do
        v := ChangeRing(e@@Ld`to_disc, Rationals());
        Append(~vs, Vector(Integers(), [Integers()!c : c in Eltseq(v)]));
    end for;
    iso := []; supp := AssociativeArray(); Qeta := AssociativeArray();
    for i->v in vs do
        q := (v*Qm, v);
        if q mod (2*dn^2) eq 0 then
            Append(~iso, i);
            Qeta[i] := q div (2*dn^2);
            supp[i] := Set(PrimeDivisors(LCM([Denominator(c/dn) : c in Eltseq(v)])));
        end if;
    end for;
    i0 := rep{i : i in iso | IsZero(vs[i])};
    Nprimes := Set(PrimeDivisors(N));
    classes_eta := [Set(s) : s in Sort([Sort(Setseq(S)) : S in {supp[i] : i in iso}])];
    sample := [];                                    // up to EtaPerClass cosets per support class
    for S in classes_eta do
        idx := [i : i in iso | supp[i] eq S];
        sample cat:= idx[1..Minimum(#idx, EtaPerClass)];
    end for;
    printf "  |L^v/L| = %o, isotropic %o, support classes %o, sampled %o\n", #elts, #iso,
           [<S, #[i : i in iso | supp[i] eq S]> : S in classes_eta], #sample;

    // G(c, a; eta) in Q[x]/(x^n - 1), by CRT over the prime powers of c
    gauss := function(c, a, v)
        G := P!1;
        for fac in Factorization(c) do
            p := fac[1]; k := fac[2]; pk := p^k; cc := c div pk;
            cnt := [0 : j in [1..pk]];
            for nu in CartesianPower([0..pk-1], 3) do
                nv := Vector(Integers(), [nu[1], nu[2], nu[3]]);
                qn := (nv*Qm, nv) div 2;                     // Q(nu)
                pr := (v*Qm, nv);                            // dn (eta, nu)
                assert pr mod dn eq 0;
                // nu = (c/p^k) nu_p + ... : the quadratic term picks up the cofactor, the linear term does not
                j := (a*cc*qn + pr div dn) mod pk;
                cnt[j + 1] +:= 1;
            end for;
            G := fold(G * &+[ P | cnt[j + 1]*zet(j*(n div pk)) : j in [0..pk-1] ]);
        end for;
        return G;
    end function;

    monos := {@ @};
    for f in fs do for r in Exponents(f) do Include(~monos, r); end for; end for;
    for r in monos do error if &+r ne 1, "a monomial of weight other than 1/2"; end for;

    // per coset representative gamma (c > 0): monomial -> <kexp, scale, m1, c0 poly>, meaning
    //   a_0(monomial | gamma) * e(1/8) = scale * sqrt(m1) * e(kexp/n) * c0
    classdata := function(g)
        tri := []; gds := [];
        for dd in ds do
            ad, bd, ed, gd := triang(g, dd); Append(~tri, <ad, bd, ed>); Append(~gds, gd);
        end for;
        W := LCM([t[3] : t in tri]);
        epsk := [ Integers()!(n*epsexp(gd)) : gd in gds ];
        bk := [ (M*tri[i][2]) div tri[i][3] : i in [1..#ds] ];     // e(b/(24 e)) = e((M b/e)/n)
        assert forall{i : i in [1..#ds] | (M*tri[i][2]) mod tri[i][3] eq 0};
        data := AssociativeArray(); nser := 0;
        for r in monos do
            L := &+[ r[i]*tri[i][1]*(W div tri[i][3]) : i in [1..#ds] ];
            if L gt 0 then continue; end if;
            if L eq 0 then
                c0poly := P!1;
            else
                nser +:= 1;
                c0, F, Wp := unit_constant(r, tri, W, L);
                if c0 eq 0 then continue; end if;
                if Type(F) eq FldRat then
                    c0poly := P!c0;
                else
                    cs := Eltseq(c0);
                    c0poly := &+[ P | cs[i]*zet((i - 1)*(n div Wp)) : i in [1..#cs] ];
                end if;
            end if;
            kexp := &+[ r[i]*(epsk[i] + bk[i]) : i in [1..#ds] ];
            E := &*[ Rationals() | tri[i][3]^r[i] : i in [1..#ds] ];   // prod e_d^{r_d}; need E^{-1/2}
            num := Numerator(E); den := Denominator(E);
            m1, s1 := SquarefreeFactorization(num*den);                 // E^{-1/2} = s1 sqrt(m1)/num
            data[r] := <kexp mod n, s1/num, m1, c0poly>;
        end for;
        return data, nser;
    end function;

    // rho_B(gamma)_{0,eta} * e(-1/8) = scale * sqrt(m0) * poly
    rhofac := function(g, i)
        a := g[1][1]; c := g[2][1]; d := g[2][2];
        poly := fold(zet(d*Qeta[i]*(n div c)) * gauss(c, a, vs[i]));
        m0, s := SquarefreeFactorization(2*c);                       // c^{-3/2}/sqrt 2 = s sqrt(m0)/(2 c^2)
        return poly, s/(2*c^2*D*N), m0;
    end function;

    // the class representatives
    Ng := func< g | (EulerPhi(M div g)*M*EulerPhi(g)) div (g*EulerPhi(M)) >;
    assert &+[Ng(g) : g in Divisors(M)] eq M*&*[Rationals() | 1 + 1/p : p in PrimeDivisors(M)];
    Smat := Matrix(Integers(), 2, 2, [0, -1, 1, 0]);
    rep1 := func< g | g eq 1 select Smat else Matrix(Integers(), 2, 2, [1, 0, g, 1]) >;
    rep2 := function(g)          // a second coset in the class, with d not = +-1 mod g when possible
        d0 := g - 1;
        for d in [2..g-2] do
            if GCD(d, g) eq 1 then d0 := d; break; end if;
        end for;
        a0 := Integers()!((Integers(g)!d0)^-1); if a0 eq 0 then a0 := 1; end if;
        return Matrix(Integers(), 2, 2, [a0, (a0*d0 - 1) div g, g, d0]);
    end function;

    t0 := Realtime();
    reps := [];                                  // <g, Ng, data> for g < M
    for g in [g : g in Divisors(M) | g lt M] do
        gm := rep1(g);
        data, nser := classdata(gm);
        Append(~reps, <gm, Ng(g), data>);
        if nser gt 0 then printf "  class %o: %o monomials needed a q-series\n", g, nser; end if;
    end for;
    // the identity class, for the zero coset only: a_0(f) at oo, rho = 1
    idtri := [<dd, 0, 1> : dd in ds];
    iddata := AssociativeArray();
    for r in monos do
        L := &+[ r[i]*ds[i] : i in [1..#ds] ];
        if L gt 0 then continue; end if;
        c0 := L eq 0 select Rationals()!1 else unit_constant(r, idtri, 1, L);
        if c0 ne 0 then iddata[r] := c0; end if;
    end for;
    printf "  class data for %o classes: %o s\n", #reps, Realtime(t0);

    // per form and sampled eta: the contributions, keyed by the squarefree part of the surd
    contrib := function(f, i, gm, Ngm, data)
        acc := AssociativeArray();
        rpoly, rscale, m0 := rhofac(gm, i);
        a0 := AssociativeArray();
        for r in Exponents(f) do
            if not IsDefined(data, r) then continue; end if;
            kexp, sc, m1, c0poly := Explode(data[r]);
            term := (f`coeffs[r] * sc) * fold(zet(kexp) * c0poly);
            if IsDefined(a0, m1) then a0[m1] +:= term; else a0[m1] := term; end if;
        end for;
        for m1 in Keys(a0) do
            mm, ss := SquarefreeFactorization(m0*m1);
            term := (Ngm * rscale * ss) * fold(rpoly * a0[m1]);
            if IsDefined(acc, mm) then acc[mm] +:= term; else acc[mm] := term; end if;
        end for;
        return acc;
    end function;
    reduce := function(acc)
        tot := P!0;
        for mm in Keys(acc) do tot +:= fold(acc[mm] * sqrtpoly(mm)); end for;
        return tot mod Phi;
    end function;

    results := []; byclass := [];
    for fi->f in fs do
        t1 := Realtime();
        vals := AssociativeArray();
        for i in sample do
            tot := P!0;
            for rp in reps do
                tot +:= reduce(contrib(f, i, rp[1], rp[2], rp[3]));
            end for;
            if i eq i0 then
                tot +:= &+[ Rationals() | f`coeffs[r]*iddata[r] : r in Exponents(f) | IsDefined(iddata, r) ];
            end if;
            tot := tot mod Phi;
            error if Degree(tot) gt 0, Sprintf("c_eta(0) is not rational at form %o, coset %o: degree %o", fi, i, Degree(tot));
            vals[i] := Coefficient(tot, 0);
        end for;
        arr := AssociativeArray(); cls := AssociativeArray();
        for S in classes_eta do
            idx := [i : i in sample | supp[i] eq S];
            vs_ := [vals[i] : i in idx];
            error if #Set(vs_) ne 1, Sprintf("c_eta(0) not constant on the support class %o: %o", S, vs_);
            cls[S] := vs_[1]/2;
            if #S eq 1 then arr[Representative(S)] := vs_[1]/2; end if;
        end for;
        Append(~results, arr); Append(~byclass, cls);
        printf "  form %o: %o  (%o s)\n", fi, [<S, cls[S]> : S in classes_eta], Realtime(t1);
        if Verify then
            for rp in reps do
                g := rp[1][2][1];
                if g eq 1 then continue; end if;             // class 1: the two reps differ by T^k, trivially equal
                gm2 := rep2(g);
                if gm2 eq rp[1] then continue; end if;
                data2 := classdata(gm2);
                for i in sample do
                    v1 := reduce(contrib(f, i, rp[1], rp[2], rp[3]));
                    v2 := reduce(contrib(f, i, gm2, rp[2], data2));
                    error if v1 ne v2, Sprintf("VERIFY FAILED: class %o, reps %o and %o differ at form %o, coset %o",
                                               g, Eltseq(rp[1]), Eltseq(gm2), fi, i);
                end for;
            end for;
            printf "    verify: every class agrees between two representatives\n";
        end if;
    end for;
    return results, byclass;
end function;

// ---------------------------------------------------------------------------------------------
// driver
// ---------------------------------------------------------------------------------------------

Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
Ld := ShimuraCurveLattice(D, N);
M := IsOdd(D*N) select 4*D*N else 2*D*N;
if lll then
    // the four weight-1/2 eta quotient MONOMIALS with the input character of composite/compmult5.m
    // (poles at every cusp, so the q-series path of every class is exercised)
    ds := Divisors(M); nd := #ds; ps := PrimeDivisors(M);
    nv := nd + 2 + #ps;
    rows := [];
    Append(~rows, [1 : d in ds] cat [0 : i in [1..2 + #ps]]);
    Append(~rows, ds cat [-24, 0] cat [0 : p in ps]);
    Append(~rows, [M div d : d in ds] cat [0, -24] cat [0 : p in ps]);
    for i->p in ps do
        Append(~rows, [Valuation(d, p) : d in ds] cat [0, 0] cat [j eq i select -2 else 0 : j in [1..#ps]]);
    end for;
    A := Matrix(Integers(), rows);
    rhs := Vector(Integers(), [1, 0, 0] cat [p eq 2 select 1 else 0 : p in ps]);
    ok, x0 := IsConsistent(Transpose(A), rhs); assert ok;
    K := KernelMatrix(Transpose(A));
    Lk := LatticeWithBasis(K);
    cl := Eltseq(ClosestVectors(Lk, -Vector(Rationals(), Eltseq(x0)) : Max := 1)[1]);
    v0 := Eltseq(x0);
    base := [Integers() | v0[i] + cl[i] : i in [1..nv]];
    xs := [base];
    B := LLL(K);
    for j in [1..3] do
        tv := Eltseq(B[j]);
        Append(~xs, [Integers() | base[i] + tv[i] : i in [1..nv]]);
    end for;
    rs := [x[1..nd] : x in xs];
    R := EtaQuotientsRing(M, D*N);
    fs := [EtaQuotient(R, r) : r in rs];
    keys := [1..#fs];
    for r in rs do printf "  LLL form r = %o\n", r; end for;
else
    fsa := BorcherdsForms(star, curves : Prec := 100);
    keys := Sort([k : k in Keys(fsa)]);
    fs := [fsa[k] : k in keys];
end if;
printf "BASE %o %o  M = %o  forms %o\n", D, N, M, keys;

if rhocheck then
    // the rho formula against the FFT, numerically, on a sample of cosets including non-isotropic ones
    CC := ComplexField(30); ii := CC.1; ee := func< z | Exp(2*Pi(CC)*ii*z) >;
    fftdata := VVWeilFFT(Ld, CC : Dual := true); elts := fftdata[7];
    Qr := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
    vs := [ChangeRing(g@@Ld`to_disc, Rationals()) : g in elts];
    nG := #elts;
    idx := [1 + (k*97) mod nG : k in [0..23]];
    mats := [ [0,-1,1,0], [1,0,2,1], [1,0,3,1], [3,1,5,2], [1,0,6,1], [7,2,10,3], [5,2,12,5], [2,1,15,8], [3,1,20,7], [1,2,4,9] ];
    for mm in mats do
        g := Matrix(Integers(), 2, 2, mm);
        if g[2][1] gt 20 then continue; end if;
        v := VVRhoInvE0FFT(fftdata, VVSTWord(g));
        a := g[1][1]; c := g[2][1]; d := g[2][2];
        // most components vanish (the Gauss sum is zero off a subgroup), so compare by DIFFERENCE with
        // the lift sign s = +-1 fixed on the first nonzero component -- a ratio of two zeros is noise
        preds := []; preds_a := [];
        for i in idx do
            q := (vs[i]*Qr, vs[i])/(2*dn^2);
            G := CC!0;
            for nu in CartesianPower([0..c-1], 3) do
                nv := ChangeRing(Vector(Integers(), [nu[1], nu[2], nu[3]]), Rationals());
                G +:= ee((a*(nv*Qr, nv)/2 + (vs[i]*Qr, nv)/dn)/c);
            end for;
            pred := ee(1/8) * CC!c^(-3/2) / (D*N*Sqrt(CC!2)) * G;
            Append(~preds, pred*ee(d*q/c));
            Append(~preds_a, pred*ee(a*q/c));
        end for;
        nzk := [k : k in [1..#idx] | Abs(preds[k]) gt 10^-10];
        if #nzk eq 0 then
            printf "RHO %o: every sampled component vanishes; max |fft| = %o\n", mm, RealField(5)!Maximum([Abs(v[i]) : i in idx]);
            continue;
        end if;
        k := nzk[1];
        s := Round(Re(v[idx[k]]/preds[k]));
        err := Maximum([Abs(v[idx[k]] - s*preds[k]) : k in [1..#idx]]);
        err_a := Maximum([Abs(v[idx[k]] - s*preds_a[k]) : k in [1..#idx]]);
        nz := #[k : k in [1..#idx] | Abs(preds[k]) gt 10^-10];
        printf "RHO %o: lift sign %o, max |fft - formula| = %o over %o cosets (%o nonzero); with a in place of d: %o\n",
               mm, s, RealField(5)!err, #idx, nz, RealField(5)!err_a;
    end for;
end if;

t0 := Realtime();
mults, byclass := exact_m0(fs, Ld, D, N : Verify := verify);
printf "EXACT multipliers: %o  (%o s)\n", [[<p, mults[i][p]> : p in Sort(Setseq(Keys(mults[i])))] : i in [1..#fs]], Realtime(t0);
if cmp then
    t0 := Realtime();
    num := M0MultipliersBySupport(fs, Ld, D, N);
    printf "M0MultipliersBySupport: %o  (%o s)\n", [[<p, num[i][p]> : p in Sort(Setseq(Keys(num[i])))] : i in [1..#fs]], Realtime(t0);
    agree := forall{i : i in [1..#fs] | forall{p : p in Keys(num[i]) | num[i][p] eq mults[i][p]}};
    printf "AGREEMENT exact vs numeric on all %o forms: %o\n", #fs, agree;
end if;
quit;
