// kappa-minus, step 2: independent local densities at the LEVEL PRIME.
//
// Convention pinned by the p=2 probe: the code's Q(x) is (1/2) x^T A x, not
// x^T A x (on the hyperbolic A = [[0,1],[1,0]] the code gives W(1) = 1/2 at odd
// m, which is the density of x1*x2 = m, not of 2*x1*x2 = m).
//
// alpha_p(m) = lim_k p^{-k} #{ x in (Z/p^k)^2 : (1/2) x^T A x = m mod p^k },
// convergence in k checked explicitly.  Compared against LocalWhittakerAtOne.
//
// The case that has never been checked independently: A = p * (unimodular),
// i.e. ord_p(det) = 2 -- the shape of L_- at a FIRING discriminant, where the
// closed form says W_p(1) = (p-1)*ord_p(m).  Both the isotropic (split p, the
// relevant one: firing forces N split) and anisotropic scaled blocks are run.
AttachSpec("ShimuraQuotients.spec");

alpha := function(p, A, m, kmax)
    prev := Rationals()!0;
    for k in [2..kmax] do
        pk := p^k; cnt := 0;
        for x1 in [0..pk-1] do
            for x2 in [0..pk-1] do
                num := A[1][1]*x1^2 + 2*A[1][2]*x1*x2 + A[2][2]*x2^2;
                if IsEven(num) and ((num div 2) - m) mod pk eq 0 then cnt +:= 1; end if;
            end for;
        end for;
        cur := Rationals()!cnt / pk;
        if k ge 3 and cur eq prev then return cur, true; end if;
        prev := cur;
    end for;
    return prev, false;
end function;

L := RSpaceWithBasis(IdentityMatrix(Integers(), 2));
mu0 := Vector([Rationals() | 0, 0]);

printf "=== calibration: p = 2, the blocks tests/Whittaker2.m validates ===\n";
for vv in [0,1] do
    A := Matrix(Integers(), 2, 2, [0, 2^vv, 2^vv, 0]);
    bad := 0; n := 0;
    for m in [1..12] do
        w := LocalWhittakerAtOne(Rationals()!m, 2, mu0, L, A);
        a, ok := alpha(2, A, m, 9);
        n +:= 1;
        if not ok or w ne a then bad +:= 1;
            printf "   v=%o m=%-3o W=%-6o alpha=%-6o %o\n", vv, m, w, a,
                ok select "MISMATCH" else "(not converged)";
        end if;
    end for;
    printf "  v = %o : %o of %o agree EXACTLY (W(1) = alpha, no constant)\n", vv, n-bad, n;
end for;

printf "\n=== the untested regime: ODD p, and p-SCALED (ord_p(det) = 2) ===\n";
// A = scale * block; blocks: hyperbolic (isotropic / p split) and diag (anisotropic)
tests := [];
for p in [3, 5, 7] do
    u := 1; while KroneckerSymbol(u, p) ne -1 do u +:= 1; end while;
    for sc in [1, p] do
        Append(~tests, <p, sc, "hyperbolic (p split)",
                        Matrix(Integers(),2,2,[0, sc, sc, 0])>);
        Append(~tests, <p, sc, "diagonal   (p inert)",
                        Matrix(Integers(),2,2,[2*sc, 0, 0, 2*sc*u])>);
    end for;
end for;

for t in tests do
    p := t[1]; sc := t[2]; nm := t[3]; A := t[4];
    kmax := p eq 3 select 6 else (p eq 5 select 4 else 3);
    mlist := [1..(p eq 7 select 10 else 14)];
    bad := 0; n := 0; noconv := 0; rows := [];
    for m in mlist do
        w := LocalWhittakerAtOne(Rationals()!m, p, mu0, L, A);
        a, ok := alpha(p, A, m, kmax);
        n +:= 1;
        if not ok then noconv +:= 1; continue; end if;
        if w ne a then bad +:= 1;
            if #rows lt 5 then Append(~rows, <m, w, a>); end if;
        end if;
    end for;
    printf "  p=%o scale=%-2o %-22o : %o of %o agree%o%o\n",
        p, sc, nm, n-bad-noconv, n-noconv,
        noconv gt 0 select Sprintf("  (%o not converged)", noconv) else "",
        bad eq 0 select "   OK" else Sprintf("   MISMATCH <m,W,alpha> = %o", rows);
end for;
printf "DONE\n";
quit;
