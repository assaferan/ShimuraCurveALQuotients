// ROUTE B: predict the Borcherds-obstruction deficit from the WEIL REPRESENTATION alone.
//
//   magma -b weildim.m
//
// Formula (Borcherds/Skoruppa shape), with alpha(A) = sum of the eigenvalue phases of A:
//
//   dim M_k(rho) - dim S_{2-k}(rho*) = d + d*k/12
//                                      - alpha( e(k/4)  * rho(S)  )
//                                      - alpha( e(-k/6) * rho(ST) )
//                                      - alpha( rho(T) )
//
// where e(x) = exp(2 pi i x), and for A of finite order with eigenvalues e(lambda_j),
// lambda_j in [0,1), alpha(A) = sum_j lambda_j.
//
// The constants were FIXED BY MATCHING THE CLASSICAL DIMENSIONS of M_k(SL_2(Z)) for the
// trivial representation (d = 1, rho(S) = rho(T) = 1), which are independent ground truth --
// NOT by tuning against this project's data.  An earlier attempt using TRACES instead of
// eigenvalue sums failed this self-test (gave 0.907 where the answer is 1), which is why the
// self-test runs first and the script refuses to report per-base numbers if it fails.
//
// alpha(A) is computed from TRACES OF POWERS rather than numerical eigenvalues:
// if A^n = I then the multiplicity of the eigenvalue e(j/n) is
//     m_j = (1/n) * sum_{m=0}^{n-1} tr(A^m) * e(-jm/n),
// which is numerically stable and integral, so it can be rounded exactly.
AttachSpec("ShimuraQuotients.spec");

PREC := 80;
CC := ComplexField(PREC);
PI := Pi(CC);
e_ := func<x | Exp(2*PI*CC.1*x)>;

// alpha(A) for a complex matrix A of finite order (order searched up to maxord).
function alphaOf(A, maxord)
    d := Nrows(A);
    Id := IdentityMatrix(CC, d);
    n := 0;
    P := A;
    for m in [1..maxord] do
        if Maximum([Abs(P[i][j] - Id[i][j]) : i in [1..d], j in [1..d]]) lt 1e-20 then
            n := m; break;
        end if;
        P := P*A;
    end for;
    error if n eq 0, Sprintf("alphaOf: matrix has no finite order <= %o", maxord);
    // traces of powers
    trs := [CC | ];
    P := IdentityMatrix(CC, d);
    for m in [0..n-1] do
        Append(~trs, Trace(P));
        P := P*A;
    end for;
    a := 0;
    for j in [0..n-1] do
        mj := (1/n) * &+[CC | trs[m+1] * e_(-j*m/n) : m in [0..n-1]];
        mjr := Round(Real(mj));
        error if Abs(mj - mjr) gt 1e-10,
              Sprintf("alphaOf: non-integral multiplicity %o at j = %o", mj, j);
        a +:= mjr * (j/n);
    end for;
    return a;
end function;

function dimMk(d, k, Sc, STc, Tc, maxord)
    a1 := alphaOf(e_(k/4)  * Sc,  maxord);
    a2 := alphaOf(e_(-k/6) * STc, maxord);
    a3 := alphaOf(Tc, maxord);
    return d + d*k/12 - a1 - a2 - a3, a1, a2, a3;
end function;

// ---------- SELF-TEST ----------
printf "=== SELF-TEST: trivial rep (d=1) vs classical dim M_k(SL_2(Z)) ===\n";
one := IdentityMatrix(CC, 1);
known := [<4,1>,<6,1>,<8,1>,<10,1>,<12,2>,<14,1>,<16,2>,<18,2>,<20,2>,<22,2>,<26,2>];
nfail := 0;
for p in known do
    k, want := Explode(p);
    got := dimMk(1, k, one, one, one, 48);
    ok := Abs(got - want) lt 1e-9;
    if not ok then nfail +:= 1; end if;
    printf "  k = %o formula = %o classical = %o  %o\n",
           k, RealField(10)!Real(got), want, ok select "ok" else "MISMATCH";
end for;
printf "SELFTEST %o (%o mismatches)\n", nfail eq 0 select "PASS" else "FAIL", nfail;
if nfail ne 0 then
    printf "ABORTING: formula is wrong, per-base numbers would be meaningless.\nDONE\n";
    quit;
end if;

// ---------- per-base ----------
// Route A ground truth: 38_5 -> deficit 1 ; 38_7 -> 0 ; 34_3 -> 0.
bases := [<34,3,0>, <38,7,0>, <38,5,1>];
printf "\n=== PER-BASE, k = 3/2 (route A truth in parentheses) ===\n";
for b in bases do
    D := b[1]; N := b[2]; truth := b[3];
    t0 := Realtime();
    Ld := ShimuraCurveLattice(D, N);
    S, T, elts, K := WeilRepresentationST(Ld);
    t_rep := Realtime() - t0;
    n0 := CyclotomicOrder(K);
    d := Nrows(S);
    printf "WEILDIM %o_%o d %o  (rep built in %os)\n", D, N, d, RealField(6)!t_rep;
    zeta := e_(1/n0);
    embm := func<M | Matrix(CC, Nrows(M), Ncols(M),
              [ &+[CC | CC!(Eltseq(M[i][j])[l]) * zeta^(l-1) : l in [1..#Eltseq(M[i][j])]]
                : i in [1..Nrows(M)], j in [1..Ncols(M)] ])>;
    Sc := embm(S); Tc := embm(T);
    STc := Sc*Tc;
    k := 3/2;
    dm, a1, a2, a3 := dimMk(d, k, Sc, STc, Tc, 48);
    printf "WEILDIM %o_%o k %o dimM %o  a1 %o a2 %o a3 %o  (truth %o) [%os]\n",
           D, N, k, RealField(10)!Real(dm), RealField(8)!a1, RealField(8)!a2,
           RealField(8)!a3, truth, RealField(6)!(Realtime() - t0);
end for;
printf "DONE\n";
quit;
