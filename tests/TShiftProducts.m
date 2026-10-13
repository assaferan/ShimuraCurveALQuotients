// tests/TShiftProducts.m
//
// The solver script's shift fallback for a tall rung with poles at infinity only: with the
// enumerated rung (380, 92, 0) and the level-380 weight-0 shift set in a cache, a request for
// (380, 158, 0) is answered by the products of the rung with the shifts (no Normaliz), every one of
// which is a lattice point of the enumerated (380, 158, 0), and the products span the same space of
// forms: the rank of the q-expansions at infinity equals that of the full rung.  The file carries
// the first-line mark of a derived spanning set.
//
// The file must hold exactly the points of the base rung and every product of a base point with a
// shift whose pole at infinity is at most 158 (nothing thinned below the size limit); with the size
// limit set to 500 it must hold, at each pole order, the first min(count, 5) of them (5 = 500 over
// the 88 pole orders present, at least four), with the warning line on the error stream.  A derived
// cached rung is never the base: with a derived file at (380, 120, 0) also in the cache, the
// request is still built from the enumerated rung 92.  NMZ_NO_SHIFT=1 bypasses the fallback.
//
// The two rungs are the committed cache files (Normaliz enumerations).  The rank is computed here
// twice with plain power series, independently of the eta-quotient code, and the two numbers are
// compared with each other, not with a stored constant.

// evaluated at top level: `eval` inside a procedure segfaults under the harness
M := 380; nlo := 92; nhi := 158;
full := Sort(eval Read(Sprintf("polymake/polymake_solution_%o_%o_0", M, nhi)));
tmp := Pipe("mktemp -d", "");
tmp := tmp[1..#tmp-1];
assert #tmp gt 0 and tmp[1] eq "/";
ok := System(Sprintf("mkdir -p %o/polymake && cp polymake/polymake_solution_%o_%o_0 polymake/tshift_w0_%o.txt %o/polymake/",
                     tmp, M, nlo, M, tmp));
assert ok eq 0;
here := GetCurrentDirectory();
rc := System(Sprintf("cd %o && NORMALIZ_BIN=/nonexistent python3 %o/nmzsolve.py %o %o 0 polymake/polymake_solution_%o_%o_0 2> nmzsolve.err",
                     tmp, here, M, nhi, M, nhi));
assert rc eq 0;                                              // served from the cache, no Normaliz
out := Sprintf("%o/polymake/polymake_solution_%o_%o_0", tmp, M, nhi);
first := Pipe(Sprintf("head -1 %o", out), "");
shifted := Sort(eval Read(out));
base := Sort(eval Read(Sprintf("polymake/polymake_solution_%o_%o_0", M, nlo)));
shifts := [t : t in eval Read(Sprintf("polymake/tshift_w0_%o.txt", M)) | exists{x : x in t | x ne 0}];
// the same request with the size limit at 500 points
rc := System(Sprintf("cd %o && rm polymake/polymake_solution_%o_%o_0 && NMZ_SHIFT_KEEP=500 NORMALIZ_BIN=/nonexistent python3 %o/nmzsolve.py %o %o 0 polymake/polymake_solution_%o_%o_0 2> capped.err",
                     tmp, M, nhi, here, M, nhi, M, nhi));
assert rc eq 0;
capped := Sort(eval Read(out));
capped_err := Read(Sprintf("%o/capped.err", tmp));
// a derived file at a rung between the two must not be used as the base
rc := System(Sprintf("cd %o && rm polymake/polymake_solution_%o_%o_0 && (echo '// derived: test'; cat polymake/polymake_solution_%o_%o_0) > polymake/polymake_solution_%o_120_0 && NORMALIZ_BIN=/nonexistent python3 %o/nmzsolve.py %o %o 0 polymake/polymake_solution_%o_%o_0 2> again.err",
                     tmp, M, nhi, M, nlo, M, here, M, nhi, M, nhi));
assert rc eq 0;
again := Sort(eval Read(out));
again_err := Read(Sprintf("%o/again.err", tmp));
// and the bypass: no fallback, so no file and a nonzero exit (Normaliz is pointed at nothing)
rc := System(Sprintf("cd %o && rm polymake/polymake_solution_%o_%o_0 && NMZ_NO_SHIFT=1 NORMALIZ_BIN=/nonexistent python3 %o/nmzsolve.py %o %o 0 polymake/polymake_solution_%o_%o_0 2> bypass.err",
                     tmp, M, nhi, here, M, nhi, M, nhi));
bypass_rc := rc; bypass_exists := FileExists(out);
_ := System(Sprintf("rm -rf %o", tmp));

procedure test_TShiftProducts(M, nlo, nhi, full, shifted, first, base, shifts, capped, capped_err, again, again_err, bypass_rc, bypass_exists)
    printf "Testing the shift fallback, (380, 92, 0) x shifts -> (380, 158, 0)...";
    ds := Divisors(M);
    assert first[1..10] eq "// derived";                         // marked as a spanning set
    assert #shifted gt 0 and #shifted lt #full;
    assert forall{r : r in shifted | r in full};                  // every product is a lattice point
    // every product with pole at most nhi is kept, and nothing else: the file is exactly the base
    // rung together with those products (one layer of shifts reaches 158 from 92: the shifts have
    // poles up to 93)
    pole := func<r | -(&+[ds[i]*r[i] : i in [1..#ds]]) div 24>;
    products := {[r[i] + t[i] : i in [1..#ds]] : r in base, t in shifts};
    expected := {r : r in base} join {v : v in products | pole(v) le nhi};
    assert Set(shifted) eq expected;
    // the size limit: at each pole order the first min(count, per) points, per = max(4, 500 div #orders)
    byorder := AssociativeArray();
    for v in Sort(SetToSequence(expected)) do
        o := pole(v);
        if not IsDefined(byorder, o) then byorder[o] := 0; end if;
        byorder[o] +:= 1;
    end for;
    per := Maximum(4, 500 div #Keys(byorder));
    assert #capped eq &+[Minimum(byorder[o], per) : o in Keys(byorder)];
    assert Set(capped) subset expected;
    assert Sprintf("keeping %o per pole order", per) in capped_err;
    // the derived file at the rung 120 is not the base
    assert "base rung 92" in again_err and Set(again) eq expected;
    // the bypass runs no fallback
    assert bypass_rc ne 0 and not bypass_exists;
    // rank of the q-expansions at infinity over GF(p): prod_d eta(d tau)^{r_d} = q^{(sum d r_d)/24} prod (1 - q^{dn})^{r_d};
    // a weight-1/2 form holomorphic away from infinity with a pole of order at most nhi there has
    // at most index/24 + nhi zeros at infinity, so this many coefficients determine it
    index := Integers()!(M * &*[1 + 1/p : p in PrimeDivisors(M)]);
    T := index div 24 + nhi + 2;
    function rank_of(pts)
        F := GF(1000003);
        R<q> := PowerSeriesRing(F, T + nhi + 1);
        E := [&*[R | 1 - q^(d*n) : n in [1..(T + nhi) div d]] : d in ds];
        rows := [];
        for r in pts do
            v24 := &+[ds[i]*r[i] : i in [1..#ds]];
            assert v24 mod 24 eq 0 and v24 div 24 ge -nhi;
            f := &*[R | E[i]^r[i] : i in [1..#ds] | r[i] ne 0];
            f := q^(v24 div 24 + nhi) * f;
            Append(~rows, [Coefficient(f, j) : j in [0..T + nhi]]);
        end for;
        return Rank(Matrix(F, rows));
    end function;
    rk_full := rank_of(full);
    rk_shift := rank_of(shifted);
    assert rk_full eq rk_shift;
    printf "Done! (%o products of %o points, rank %o = %o; %o kept under the size limit)\n", #shifted, #full, rk_shift, rk_full, #capped;
end procedure;

test_TShiftProducts(M, nlo, nhi, full, shifted, first, base, shifts, capped, capped_err, again, again_err, bypass_rc, bypass_exists);
