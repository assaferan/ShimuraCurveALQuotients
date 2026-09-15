// Regenerates the expected-value table that tests/Whittaker2.m parts (4) and (5) commit.
//
// WHY IT IS OFFLINE. It brute-forces a 2-adic representation density by literally enumerating
// (Z/2^k)^2, which is 2^(2k) Magma loop iterations per configuration -- about 12 s at k = 12 on this
// Mac, and 16x that for each further rung. tests/Whittaker2.m must stay in the 1 s ballpark because
// the suite runs it every time, so it commits the CONVERGED VALUES and this file, hand-run, is what
// re-derives them. Files in tests/_offline/ are invisible to run_tests.m's `ls tests/*.m` sweep and
// to the CI matrix (both filter a leading `_`).
//
//     W2O_KMAX=12 magma -b filename:=tests/_offline/Whittaker2Oracle.m run_tests.m < /dev/null > out.txt
//
// ⚠ ALWAYS REDIRECT STDIN; and Magma buffers stdout to a file, so a log that is not growing is not
// evidence that nothing is happening. Full sweep at the default kmax = 12: 25 Gram matrices,
// 136 configurations, roughly 30-35 min.
//
// WHAT IT COMPUTES, and the convention that matters.  In this repo Q(x) = 1/2 * x G x^T, so
//
//     alpha_2(m, mu + L) = lim_{k->oo} 2^(-k) * #{ x in mu + (Z/2^k)^2 : Q(x) = m mod 2^k },
//
// and with mu = (P/2^e, R/2^e), x = (P + 2^e y1, R + 2^e y2)/2^e, the condition Q(x) = m mod 2^k is
//     w G w^T = 2*m*4^e   mod 2^(k+1+2e),     w = (P + 2^e y1, R + 2^e y2),
// counted over (y1,y2) in (Z/2^k)^2 and divided by 2^k. Sanity check of the convention:
// G = [[0,1],[1,0]] must give Q(x,y) = xy, the standard hyperbolic plane.
//
// This is independent of BOTH Yang's explicit formula AND the library -- that is the point. The
// numbers it prints are compared against LocalWhittakerAtOne here, and printed in the exact source
// format of tests/Whittaker2.m's `production_shapes` so the committed table can be diffed against a
// fresh run.
//
// *** THE CONVERGENCE TRAP. *** The approximants PLATEAU AT A WRONG VALUE before converging, so
// never stop on "the last two agree". Measured at 6_1, d = -24 (shape d1d2), k = 1..14:
//     m = 8    1 1 1 1 1 2 2 2 2 2 2 2 2 2        m = 32   1 1 1 1 1 1 1 2 2 2 2 2 2 2
// A two-repeats rule returns 1 for both; the true value is 2 in both cases. Over the 300 mu = 0
// comparisons the LAST k at which any approximant still moved was k = 8 (re-confirmed to k = 16),
// and over the 666 nonzero-coset comparisons it was k = 2. Hence the default kmax = 12, and hence
// the whole k-sequence is printed: read the plateau, do not trust a stopping rule.

// kmax. Default 12: the last k at which ANY approximant below still moved was 8 (mu = 0) and 2
// (nonzero cosets), so 12 clears the plateau with four spare rungs; the committed table in
// tests/Whittaker2.m was taken at k = 14 and re-confirmed at k = 16.
//
// ⚠ Set it with the ENVIRONMENT, `W2O_KMAX=10 magma -b filename:=... run_tests.m < /dev/null`.
// A command-line `kmax:=10` does NOT reach here: run_tests.m loads a test with Read() and runs it
// through eval, and Magma's eval scope does not see the top-level command-line assignments. The
// variable is simply unassigned and the default silently wins -- which looks exactly like the flag
// having been honoured, only slower.
w2o_kmax := (GetEnv("W2O_KMAX") ne "") select StringToInteger(GetEnv("W2O_KMAX")) else 12;

// alpha_k(m) for k = 1..kmax, for every m at once (one histogram pass per k).
function w2o_brute(G, e, P, R, ms, kmax)
    pe := 2^e;
    out := [[Rationals() | ] : m in ms];
    for k in [1..kmax] do
        M := 2^(k+1+2*e); n := 2^k;
        Rm := Integers(M);
        a := Rm!G[1,1]; b := Rm!(2*G[1,2]); c := Rm!G[2,2];
        w1 := [Rm!(P + pe*y) : y in [0..n-1]];
        w2 := [Rm!(R + pe*y) : y in [0..n-1]];
        w2sq := [x^2 : x in w2];
        hist := [0 : i in [1..M]];
        for i in [1..n] do
            t1 := a*w1[i]^2; bw := b*w1[i];
            for j in [1..n] do
                v := Integers()!(t1 + bw*w2[j] + c*w2sq[j]);
                hist[v+1] +:= 1;
            end for;
        end for;
        for t in [1..#ms] do
            tgt := 2*ms[t]*pe^2;
            error if Denominator(tgt) ne 1, "target 2*m*4^e is not integral -- e too small for this mu";
            Append(~out[t], hist[(Integers()!tgt mod M) + 1] / n);
        end for;
    end for;
    return out;
end function;

// Kept in step with tests/Whittaker2.m's copy; see the comment there for why even-H and even-A must
// not share a label.
function w2o_shape(G)
    v := func< x | x eq 0 select 10^6 else Valuation(x, 2) >;
    s  := Minimum([v(G[1,1]), v(G[1,2]), v(G[2,2])]);
    dg := Minimum([v(G[1,1]), v(G[2,2])]);
    if dg gt s then
        u := Determinant(G) div 4^s;
        error if not (u mod 8 in {3, 7}), "not an even binary 2-adic lattice";
        return Sprintf("even%o%o", s, (u mod 8 eq 7) select "H" else "A");
    end if;
    return Sprintf("d%od%o", s, Valuation(Determinant(G), 2) - s);
end function;

function w2o_cosets(G)
    t := Valuation(Determinant(G), 2); pt := 2^t;
    GQ := ChangeRing(G, Rationals());
    return [w : w in [Vector([Rationals() | a/pt, b/pt]) : a in [0..pt-1], b in [0..pt-1]]
              | &and[IsIntegral(x) : x in Eltseq(w*GQ)]];
end function;

// The same Gram matrices tests/Whittaker2.m commits, with the (base, d) each was harvested from.
w2o_grams := [
    <"6_1",  -4,  [-6,-114,-114,-2172]>,        <"14_1", -4,  [-1414,-2968,-2968,-6230]>,
    <"10_1", -20, [-326,370,370,-420]>,         <"21_1", -4,  [-42,0,0,-42]>,
    <"26_1", -20, [-78,416,416,-2262]>,         <"34_1", -20, [-33558,10574,10574,-3332]>,
    <"38_1", -4,  [-52630,-92454,-92454,-162412]>, <"38_1", -20, [-342,-2470,-2470,-17860]>,
    <"10_3", -20, [-486,-2310,-2310,-10980]>,   <"6_5",  -4,  [-870,-510,-510,-300]>,
    <"6_1",  -24, [-22,-52,-52,-124]>,          <"14_1", -8,  [-714,-1568,-1568,-3444]>,
    <"10_1", -8,  [-5930,2600,2600,-1140]>,     <"26_1", -8,  [-26,234,234,-2158]>,
    <"26_1", -24, [-650,988,988,-1508]>,        <"34_1", -24, [-65450,2108,2108,-68]>,
    <"38_1", -24, [-266,-836,-836,-2660]>,      <"10_3", -8,  [-9690,-15000,-15000,-23220]>,
    <"6_5",  -24, [-10,-220,-220,-4900]>,
    // even1A, the 2-modular NON-hyperbolic block 2*[[2,1],[1,2]] -- what the ODD discriminants give
    <"6_1",  -3,  [-964,-62,-62,-4]>,           <"14_1", -11, [-252,-1694,-1694,-11396]>,
    <"10_1", -3,  [-9780,1170,1170,-140]>,      <"26_1", -19, [-52,1430,1430,-39572]>,
    <"34_1", -11, [-139740,5338,5338,-204]>,    <"38_1", -19, [-20,-2,-2,-4]>
];

procedure w2o_run(kmax)
    L := RSpaceWithBasis(IdentityMatrix(Integers(), 2));
    ms0 := [1,2,3,4,5,6,7,8,12,16,24,32];

    // CALIBRATE FIRST, always: reproduce the values tests/Whittaker2.m parts (1)-(3) already tie to
    // Yang Thm 4.1 before believing anything the counter says about a shape nobody has checked.
    printf "CALIBRATION on the two hyperbolic blocks (must be exact, ratio 1):\n";
    ncal := 0;
    for cal in [ <Matrix(Integers(),2,2,[0,1,1,0]), 0, 0, 0>,
                 <Matrix(Integers(),2,2,[0,2,2,0]), 0, 0, 0>,
                 <Matrix(Integers(),2,2,[0,2,2,0]), 1, 1, 1>,     // mu = (1/2,1/2)
                 <Matrix(Integers(),2,2,[0,2,2,0]), 1, 1, 0> ] do // mu = (1/2,0)
        Gc, e, P, R := Explode(cal);
        mu := Vector([Rationals() | P/2^e, R/2^e]);
        cms := [Rationals() | 1,2,3,4,5,6,8,12];
        seqs := w2o_brute(Gc, e, P, R, cms, Minimum(kmax, 12));
        for t in [1..#cms] do
            lib := LocalWhittakerAtOne(cms[t], 2, mu, L, Gc);
            brute := seqs[t][#seqs[t]];
            printf "  G=%o mu=%o m=%o: lib %o brute %o %o\n", Eltseq(Gc), Eltseq(mu), cms[t],
                   lib, brute, (lib eq brute) select "OK" else "*** MISMATCH ***";
            error if lib ne brute, "calibration failed -- the counter is wrong, fix it before going on";
            ncal +:= 1;
        end for;
    end for;
    printf "CALIBRATION: %o/%o exact.\n\n", ncal, ncal;

    n := 0; bad := 0;
    for row in w2o_grams do
        base, d, gent := Explode(row);
        G := MatrixAlgebra(Integers(), 2) ! gent;
        shape := w2o_shape(G);
        t0 := Cputime();

        seqs := w2o_brute(G, 0, 0, 0, [Rationals()!m : m in ms0], kmax);
        w0 := [];
        for t in [1..#ms0] do
            printf "%o d=%o %o mu=0 m=%o:", base, d, shape, ms0[t];
            for x in seqs[t] do printf " %o", x; end for; printf "\n";
            val := seqs[t][#seqs[t]];
            Append(~w0, val);
            lib := LocalWhittakerAtOne(Rationals()!ms0[t], 2, Vector([Rationals()|0,0]), L, G);
            n +:= 1; if lib ne val then bad +:= 1; printf "  *** MISMATCH: lib %o\n", lib; end if;
        end for;

        rows := [];
        for mu in w2o_cosets(G) do
            if IsZero(mu) then continue; end if;
            den := LCM([Denominator(x) : x in Eltseq(mu)]);
            e := Valuation(den, 2); error if 2^e ne den, "coset denominator is not a power of 2";
            P := Integers()!(mu[1]*2^e); R := Integers()!(mu[2]*2^e);
            Qmu := (Matrix(mu)*ChangeRing(G, Rationals())*Transpose(Matrix(mu)))[1,1] / 2;
            frac := Qmu - Floor(Qmu);
            error if frac eq 0, "Q(mu) is 2-integral -- the m list below would be the integral one";
            ms := [frac + j : j in [0,1,2,3,4,8]];
            mseqs := w2o_brute(G, e, P, R, ms, kmax);
            vals := [];
            for t in [1..#ms] do
                printf "%o d=%o %o mu=%o m=%o:", base, d, shape, Eltseq(mu), ms[t];
                for x in mseqs[t] do printf " %o", x; end for; printf "\n";
                val := mseqs[t][#mseqs[t]];
                Append(~vals, val);
                lib := LocalWhittakerAtOne(ms[t], 2, mu, L, G);
                n +:= 1; if lib ne val then bad +:= 1; printf "  *** MISMATCH: lib %o\n", lib; end if;
            end for;
            Append(~rows, Sprintf("<[Rationals()|%o,%o], %o>", mu[1], mu[2], vals));
        end for;

        // emitted in tests/Whittaker2.m's `production_shapes` source format, for a direct diff
        printf "\n  < \"%o\", %o, \"%o\", %o,\n    %o,\n    [ %o ] >,\n\n",
               base, d, shape, gent, w0, Join(rows, ", ");
        printf "// ^ %o d=%o done in %o s\n\n", base, d, Cputime(t0);
    end for;
    printf "\nTOTAL: %o comparisons against LocalWhittakerAtOne, %o mismatches (kmax = %o).\n",
           n, bad, kmax;
    error if bad ne 0, "brute-force oracle disagrees with the library -- do NOT edit library code; report it";
end procedure;

w2o_run(w2o_kmax);
