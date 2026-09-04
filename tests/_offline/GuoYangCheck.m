// Shared checker for Guo-Yang CM-value tables, used by the per-base tests in this directory.
//
// SOURCE. J.-W. Guo and Y. Yang, "Equations of hyperelliptic Shimura curves", arXiv:1510.06193,
// appendix tables "CM-values of X_0^D(N)" -- the primary hauptmodule column (full star quotient
// X_0(D,N)/W_{D,N}). Extracted mechanically from the arXiv v1 TeX source; see
// vvdata/weyl-campaign/guoyang-tables.md on the campaign branch for the extraction method and
// the full per-base coverage/exclusion record.
//
// WHY THESE LIVE HERE AND NOT IN tests/. Each base costs a full Borcherds + Schofer run --
// seconds for the small ones, but multiple HOURS for the heaviest (87_1 took ~2h on a 32-core
// machine). Files in tests/_offline/ are invisible to both the default `ls tests/*.m` sweep in
// run_tests.m and the CI matrix generator (.github/workflows/tests.yml filters out anything
// starting with `_`, and this whole directory does), so they never run automatically -- run one
// explicitly with `magma -b filename:=tests/_offline/GuoYang_D_N.m run_tests.m`.
//
// THE CHECK. A Hauptmodul is only defined up to a Mobius transformation, and our normalisation
// need not match Guo-Yang's. So: pick 3 discriminants the table lists, use OUR computed values
// there as one Mobius frame and Guo-Yang's PUBLISHED values there as another, and compare the
// resulting cross-ratios at every other listed discriminant -- a quantity invariant under either
// side's choice of coordinate on P^1(Q). This generalises tests/ExternalCMValues.m's single-value
// check (same cross-ratio idea, same `mobius` formula) to a whole table's rows at once.
//
// DISCRIMINANTS EXCLUDED FROM A TABLE, on principle: any discriminant Guo-Yang's table lists
// MORE THAN ONCE with different values (a real phenomenon -- two distinct CM points can share a
// discriminant -- seen at X_0^10(19), disc -52). A lookup keyed by discriminant alone can't tell
// the two points apart, so such rows are left out rather than risk pinning the wrong one.

function mobius(z0, z1, z2, z)
    // cross-ratio (z, z0; z1, z2), written so each factor drops out when its argument is Infinity
    num := (z eq Infinity() or z0 eq Infinity()) select 1 else z - z0;
    den := (z eq Infinity() or z1 eq Infinity()) select 1 else z - z1;
    c1  := (z2 eq Infinity() or z1 eq Infinity()) select 1 else z2 - z1;
    c2  := (z2 eq Infinity() or z0 eq Infinity()) select 1 else z2 - z0;
    if den*c2 eq 0 then return Infinity(); end if;
    return (num*c1)/(den*c2);
end function;

// gy: sequence of <disc, published value> pairs, value in Rationals() or Infinity().
procedure test_gy_table(D, N, gy)
    printf "  Guo-Yang CM-values of X0^%o(%o), primary hauptmodule...", D, N;
    Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
    Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
    curves := GetQuotientsAndGenera([Xstar]);
    assert exists(star){c : c in curves | IsStarCurve(c)};

    keep := { t[1] : t in gy };
    tab := ValuesAtCMPoints(star, curves : Keep := keep);
    discs := tab`Discs;
    srow := tab`Values[tab`sIndex];
    idx := AssociativeArray();
    for i->d in discs do idx[d] := i; end for;

    have := [t : t in gy | IsDefined(idx, t[1])];
    error if #have lt 4,
        Sprintf("X0^%o(%o): only %o of %o published discriminants reached our table",
                D, N, #have, #gy);

    ref := [];
    for t in have do
        if #ref eq 3 then break; end if;
        if forall{r : r in ref | r[2] ne t[2]} then Append(~ref, t); end if;
    end for;
    error if #ref lt 3, Sprintf("X0^%o(%o): could not find 3 rows with distinct published values", D, N);
    z0 := srow[idx[ref[1][1]]]; z1 := srow[idx[ref[2][1]]]; z2 := srow[idx[ref[3][1]]];
    w0 := ref[1][2]; w1 := ref[2][2]; w2 := ref[3][2];

    nchecked := 0;
    for t in have do
        d, expected := Explode(t);
        if d in {ref[1][1], ref[2][1], ref[3][1]} then continue; end if;
        got := mobius(z0, z1, z2, srow[idx[d]]);
        want := mobius(w0, w1, w2, expected);
        error if got ne want,
            Sprintf("X0^%o(%o): disc %o got %o, want %o (published %o; reference discs %o,%o,%o)",
                    D, N, d, got, want, expected, ref[1][1], ref[2][1], ref[3][1]);
        nchecked +:= 1;
    end for;
    printf " ok (%o of %o published values)\n", nchecked, #gy;
end procedure;
