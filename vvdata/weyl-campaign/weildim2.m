// ROUTE B, matrix-free.  magma -b weildim2.m
//
// WHY MATRIX-FREE: |disc group| is 20808 / 72200 / 141512 at 34_3 / 38_5 / 38_7, so the d x d
// Weil matrices (~2e10 entries at 38_7) cannot be built at all.  WeilRepresentationST is only
// usable for tiny bases.  But the dimension formula needs only TRACES OF POWERS, and every one
// of them collapses to an O(d) sum over the discriminant form:
//
//   Sum(a) := sum_gamma e(a*Q(gamma)),   G := Sum(1) (Milgram),   t2 := #{gamma : 2 gamma = 0}
//
// With S_ij = (sigma/sqrt d) e(-b(i,j)), T_ii = e(Q_i), b(i,j) = (gamma_i, gamma_j), sigma = e(1/8):
//   rho(S)^2 = sigma^2 * P   (P: gamma -> -gamma)   [character orthogonality]
//   tr S   = (sigma/sqrt d) Sum(-2)      tr S^2 = sigma^2 t2
//   tr S^3 = (sigma^3/sqrt d) Sum(2)     tr S^4 = sigma^4 d = -d          (period 8)
//   tr ST  = (sigma/sqrt d) Sum(-1)
//   tr (ST)^2 = (sigma^2/d) G Sum(-3)    [complete the square: Q(d)-2b(i,d) = Q(d-2g_i)-4Q(g_i)]
//   tr (ST)^3 = sigma^2 t2               tr (ST)^4 = (sigma^3/sqrt d) Sum(3)
//   tr (ST)^5 = (sigma^4/d) G^2          tr (ST)^6 = sigma^4 d            (period 12)
//
// The (ST)^3 = S^2 step uses (ST)^3 = -I = S^2 in SL_2(Z); the metaplectic lift is ASSUMED to
// agree, which is exactly what the 15_2 cross-check below tests.
AttachSpec("ShimuraQuotients.spec");

PREC := 40;
CC := ComplexField(PREC);
PI := Pi(CC);
e_ := func<x | Exp(2*PI*CC.1*x)>;

// alpha(A) from traces of powers: if A^n = I, multiplicity of e(j/n) is
// (1/n) sum_m tr(A^m) e(-jm/n).
function alphaFromTraces(trs, n)
    a := 0;
    for j in [0..n-1] do
        mj := (1/n) * &+[CC | trs[m+1] * e_(-j*m/n) : m in [0..n-1]];
        mjr := Round(Real(mj));
        error if Abs(mj - mjr) gt 1e-6,
              Sprintf("alphaFromTraces: non-integral multiplicity %o at j=%o", mj, j);
        a +:= mjr * (j/n);
    end for;
    return a;
end function;

// Traces of powers of rho(S) and rho(ST), and alpha(rho(T)), all in O(d).
function weilData(D, N : dual := false)
    Ld := ShimuraCurveLattice(D, N);
    dg := Ld`disc_grp; Q := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
    elts := [g : g in dg];
    d := #elts;
    vs := [ChangeRing(g@@Ld`to_disc, Rationals()) : g in elts];
    nm := [ (vs[i]*Q, vs[i])/(2*dn^2) : i in [1..d] ];
    if dual then nm := [-x : x in nm]; end if;
    sig := dual select e_(-1/8) else e_(1/8);

    Sum := func<a | &+[CC | e_(a*nm[i]) : i in [1..d]]>;
    S1 := Sum(1); Sm1 := Sum(-1); S2 := Sum(2); Sm2 := Sum(-2); S3 := Sum(3); Sm3 := Sum(-3);
    G := S1;
    zero := dg!0;
    t2 := #[i : i in [1..d] | 2*elts[i] eq zero];
    rd := Sqrt(CC!d);

    trS := [CC | d, (sig/rd)*Sm2, sig^2*t2, (sig^3/rd)*S2 ];
    trS cat:= [ sig^4*trS[i] : i in [1..4] ];                       // period 8
    trST := [CC | d, (sig/rd)*Sm1, (sig^2/d)*G*Sm3, sig^2*t2,
                  (sig^3/rd)*S3, (sig^4/d)*G^2 ];
    trST cat:= [ sig^4*trST[i] : i in [1..6] ];                     // period 12

    alphaT := &+[ Rationals() | (x - Floor(x)) where x := nm[i] : i in [1..d] ];
    return d, trS, trST, alphaT, nm, sig;
end function;

function dimMk(D, N, k : dual := false)
    d, trS, trST, alphaT := weilData(D, N : dual := dual);
    a1 := alphaFromTraces([ e_(k*m/4)  * trS[m+1]  : m in [0..7]  ], 8);
    a2 := alphaFromTraces([ e_(-k*m/6) * trST[m+1] : m in [0..11] ], 12);
    return d + d*k/12 - a1 - a2 - alphaT, d, a1, a2, alphaT;
end function;

// ---------- CROSS-CHECK the O(d) traces against explicit matrices (small base only) ----------
// 6_1 has d = 72, so the explicit matrix powers are instant.  (15_2 has d = 1800: 1800^3 per
// product, which is the very cost this route exists to avoid.)
printf "=== CROSS-CHECK O(d) formulas vs explicit Weil matrices at 6_1 (d = 72) ===\n";
D0 := 6; N0 := 1;
Ld := ShimuraCurveLattice(D0, N0);
Sm, Tm, elts0, K0 := WeilRepresentationST(Ld);
n0 := CyclotomicOrder(K0);
zeta := e_(1/n0);
embm := func<M | Matrix(CC, Nrows(M), Ncols(M),
          [ &+[CC | CC!(Eltseq(M[i][j])[l]) * zeta^(l-1) : l in [1..#Eltseq(M[i][j])]]
            : i in [1..Nrows(M)], j in [1..Ncols(M)] ])>;
Sc := embm(Sm); Tc := embm(Tm); STc := Sc*Tc;
dd, trS, trST, alphaT := weilData(D0, N0);
P := Sc; okall := true;
for m in [1..4] do
    want := Trace(P);
    got := trS[m+1];
    ok := Abs(want - got) lt 1e-8 * (1 + Abs(want));
    if not ok then okall := false; end if;
    printf "  tr rho(S)^%o  matrix %o   formula %o   %o\n", m,
           ComplexField(8)!want, ComplexField(8)!got, ok select "ok" else "MISMATCH";
    P := P*Sc;
end for;
P := STc;
for m in [1..6] do
    want := Trace(P);
    got := trST[m+1];
    ok := Abs(want - got) lt 1e-8 * (1 + Abs(want));
    if not ok then okall := false; end if;
    printf "  tr rho(ST)^%o matrix %o   formula %o   %o\n", m,
           ComplexField(8)!want, ComplexField(8)!got, ok select "ok" else "MISMATCH";
    P := P*STc;
end for;
printf "CROSSCHECK %o\n\n", okall select "PASS" else "FAIL";

// ---------- per-base ----------
printf "=== dim M_{3/2} per base (route A measured deficit in parens) ===\n";
for b in [<34,3,0>, <38,5,1>, <38,7,0>] do
    D := b[1]; N := b[2]; truth := b[3];
    for du in [false, true] do
        t0 := Realtime();
        v, d, a1, a2, aT := dimMk(D, N, 3/2 : dual := du);
        printf "WEILDIM %o_%o dual %o d %o dimM %o (a1 %o a2 %o aT %o) truth %o [%os]\n",
               D, N, du, d, RealField(10)!Real(v), RealField(8)!a1, RealField(8)!a2,
               RealField(8)!aT, truth, RealField(5)!(Realtime()-t0);
    end for;
end for;
printf "DONE\n";
quit;
