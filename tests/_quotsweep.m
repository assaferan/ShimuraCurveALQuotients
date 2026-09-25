// tests/_quotsweep.m -- NOT a test (leading underscore: excluded from the CI matrix, see
// run_tests.m).  Hand-run:
//
//     magma -b tests/_quotsweep.m < /dev/null > sweep.log
//     magma -b bases:=39_2,6_7 tests/_quotsweep.m < /dev/null
//
// WHAT IT DOES.  For every base with an X0_ test: take the top curve and the stored involutions,
// close them into the FULL Atkin-Lehner group by matrix products, quotient by each element, and
// ask which committed cover key the result matches.  It reports three things:
//   * agreements  -- the derived quotient IS the committed curve at key {1,m};
//   * mismatches  -- and, crucially, WHICH OTHER KEY the derived curve does match, because a
//     permutation of keys is a labelling problem and a match with nothing is a derivation problem;
//   * derived-with-no-expected -- the coverage this would ADD if the values were installed.
//
// MEASURED 2026-09-24 over all 50 bases:  75 matched | 4 mismatched | 112 new | 0 failed |
// 9 skipped for a non-hyperelliptic (CRV) top.  The 112 is the one-step coverage gain.
//
// ⚠⚠ ALL FOUR MISMATCHES ARE AT 39_2 AND FORM ONE TRANSPOSITION:
//        w_2  -> key {1,26}      w_26 -> key {1,2}
//        w_6  -> key {1,78}      w_78 -> key {1,6}
//    while w_3 -> {1,3} and w_13 -> {1,13} agree.  Over F_2 with basis (2,3,13) the discrepancy
//    fixes 3 and 13 and sends 2 -> 26 = 2*13, i.e. exactly w_2 <-> w_26.  Same SHAPE as the
//    settled w_10 <-> w_13 swap at 10_13.  UNRESOLVED as of 2026-09-24: either the ws_data labels
//    in tests/_offline/X0_39_2.m are wrong or the pipeline's cover-key labelling is.  The arbiter
//    is Ogg's -- the fixed points of w_m are the CM points of discriminant -4m.
//    ⚠ {1,2} and {1,26} are NON-isomorphic over Q yet share point counts at ten primes, so a
//    point-count check reads as agreement and a quadratic-twist scan finds nothing.  Only
//    IsIsomorphic separates them.
//
// ⚠ A CRV top is skipped, not failed: the recipe is written for y^2 = f(x).  Those nine bases
// (10_13 10_19 14_3 21_2 26_3 57_1 6_17 82_1 93_1) are the gap between a 76% and a 69% ceiling.

import "tests/_quotbyinvol.m" : QuotientByInvolution, QuotientMatches, ALMatrixGroupFromGenerators;

if not assigned bases then
    bases := "10_11,10_13,10_23,10_3,10_7,134_1,146_1,14_1,14_3,14_5,15_1,15_2,194_1,206_1,21_2,"
      cat "22_3,22_5,26_1,26_3,33_1,34_1,35_1,38_1,39_1,46_1,51_1,55_1,57_1,58_1,62_1,69_1,6_11,"
      cat "6_13,6_17,6_19,6_29,6_31,6_37,6_5,6_7,74_1,82_1,86_1,94_1,10_19,111_1,21_1,39_2,87_1,93_1";
end if;

// Read cover_data / ws_data out of a test WITHOUT running its body -- re-deriving a base costs
// minutes to hours, and this tool only wants the stored data.
function LoadTestData(base)
    try
        txt := Read("tests/X0_" cat base cat ".m");
    catch err
        txt := Read("tests/_offline/X0_" cat base cat ".m");
    end try;
    i := Index(txt, "end function;");
    keep := "";
    for ln in Split(txt[1..i+12], "\n") do
        // drop `import`: it compiles another package eagerly and defeats AttachSpec (CLAUDE.md)
        if #ln ge 6 and ln[1..6] eq "import" then continue; end if;
        keep cat:= ln cat "\n";
    end for;
    F := eval (keep cat "\nreturn load_covers_and_ws_data_" cat base cat ";");
    return F();
end function;

nMatch := 0; nMis := 0; nNew := 0; nFail := 0; nCRV := 0;
misrows := [];

for base in Split(bases, ",") do
    cover_data, ws_data := LoadTestData(base);
    if not IsDefined(ws_data, {1}) or #Keys(ws_data[{1}]) eq 0 then
        printf "%-8o SKIP (no ws_data at the top key)\n", base; continue;
    end if;
    top := cover_data[{1}][1];
    if Type(top) eq List then top := top[1]; end if;      // several published candidates
    if Type(top) ne CrvHyp then
        printf "%-8o SKIP (non-hyperelliptic top)\n", base; nCRV +:= 1; continue;
    end if;
    ftop := HyperellipticPolynomials(SimplifiedModel(top));
    g := (Degree(ftop) - 1) div 2;

    gens := AssociativeArray();
    for m in Keys(ws_data[{1}]) do gens[m] := ChangeRing(ws_data[{1}][m], Rationals()); end for;
    all, consistent, why := ALMatrixGroupFromGenerators(gens, g+1);
    printf "%-8o genus %-2o  %o stored -> group of %o%o\n", base, Genus(top),
           #Keys(gens), #Keys(all), consistent select "" else ("   INCONSISTENT: " cat why);

    for m in Sort(SetToSequence(Keys(all))) do
        okd, Fq, note := QuotientByInvolution(ftop, all[m]);
        if not okd then
            printf "    w_%-5o DERIVATION FAILED: %o\n", m, note; nFail +:= 1; continue;
        end if;
        key := {1, m};
        if not IsDefined(cover_data, key) then nNew +:= 1; continue; end if;
        Ck := cover_data[key][1];
        if Type(Ck) ne CrvHyp then continue; end if;
        same, why2 := QuotientMatches(Fq, Ck);
        if same then nMatch +:= 1; continue; end if;
        nMis +:= 1;
        // WHICH OBJECT is it, then?  A derived curve matching a DIFFERENT key is a labelling
        // problem; matching nothing is a derivation problem.  They need opposite responses.
        hit := "";
        for k in Keys(cover_data) do
            if k eq key then continue; end if;
            Cj := cover_data[k][1];
            if Type(Cj) ne CrvHyp then continue; end if;
            h, _ := QuotientMatches(Fq, Cj);
            if h then hit cat:= Sprintf(" %o", Sort([Integers()|q : q in k])); end if;
        end for;
        printf "    w_%-5o MISMATCH at key %o  [%o]%o\n", m,
               Sort([Integers()|q : q in key]), why2,
               hit eq "" select "  (matches NO committed key)" else ("  but MATCHES:" cat hit);
        Append(~misrows, base cat " w_" cat IntegerToString(m));
    end for;
end for;

printf "\n==== matched %o | mismatched %o | derived-with-no-expected %o | failed %o | CRV skipped %o ====\n",
       nMatch, nMis, nNew, nFail, nCRV;
if #misrows gt 0 then printf "mismatches: %o\n", misrows; end if;
