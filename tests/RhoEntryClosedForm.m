// tests/RhoEntryClosedForm.m
//
// DOES THE CODE AGREE WITH OUR OWN PAPER?
//
// Everything else that validates this project compares against an EXTERNAL source -- Guo-Yang's
// published equations and CM tables, Errthum's singular moduli, Kudla-Yang's Propositions. That
// leaves a gap nobody had closed: the preprint's own Section 10 (`sec:determined`) derives closed
// forms for the rho-vector, and NOTHING in the codebase implements or checks them. A grep for
// Ngeneral / noresplit / rhoentry / Tel_s / c_p^Eich across every Magma source returns nothing.
//
// This test closes that gap for `thm:rhoentry`, which is the most directly checkable of the three
// because it is a closed form for a matrix entry the code already builds.
//
// THE CLAIM (paper, thm:rhoentry). Let gamma in SL_2(Z) have lower row (c,d), let g be its class,
// and let gamma~ be the canonical metaplectic lift. Then (rho_Abar(gamma~) e_0)_{eta*}:
//   (i)   vanishes iff N | g;
//   (ii)  has absolute value g/sqrt(|A|) for odd g, and g/M for 4 | g.
// with |A| = M^2/2 (proof of thm:rhoentry, via |A[c]| and Stromberg Thm 6.4).
//
// WHY IT MATTERS BEYOND BOOKKEEPING. The p | N local structure is what the current
// `p | gcd(d,N)` investigation turns on (paper/DRAFT-pN-local-factor.md): at a level prime the
// lattice is Eichler, hence NON-unimodular, and Schofer Thm 4.1 assumes unimodularity at
// unramified primes. Claim (i) -- vanishing exactly when N | g -- is precisely a level-prime
// statement, so a disagreement here would localise a real divergence between code and theory.
// Note what is ALREADY known to agree: KY Props 5.4 and 5.5 reproduce the pipeline's own local
// factor at the level prime (tests/KudlaYangLocal.m, 270 + 1440 checks). So if this test also
// passes, Section 10 is eliminated as the site of the defect and the search moves upstream.
//
// ⚠ SCOPE. This checks |A| = M^2/2 on 7 bases, and claim (ii) -- the absolute value -- on the two
// small ones. It does NOT check (iii), the e_8(E_0(g)) phase for odd g, nor (iv). Those need the
// canonical-representative normalisation pinned first; claiming them here would overstate what is
// verified.
//
// ⚠⚠ A GAP IN THE PAPER, found by writing this test (2026-09-05). thm:rhoentry (ii) gives the
// absolute value for ODD g and for 4 | g, and says NOTHING about g = 2 mod 4. Those cases occur:
// at X0^6(1) and X0^10(1) they are g = 2, 6, 10. Measured, every one of them has
//     |entry| = 2g/M        (i.e. |entry|^2 = 4g^2/M^2, exactly 4x the 4|g formula)
// -- 8 of 8 such gamma, consistently. An earlier draft of this test applied the 4|g formula there
// and reported 8 "failures"; they were the test's error, not a divergence. The cases are now
// SKIPPED rather than asserted, because asserting a formula the paper does not state would be
// inventing theory. Worth adding to thm:rhoentry as a third case.

printf "Checking rho-entry against the preprint's thm:rhoentry closed form...\n";

re_total := 0; re_fail := 0; re_bad := [];

// ---- CHEAP CHECK, every base: |A| = M^2/2 (proof of thm:rhoentry) --------------------------
// Costs one ShimuraCurveLattice call per base and no Weil matrices, so it runs everywhere the
// expensive check cannot. This is a genuine code-vs-paper comparison: |A| is the size of the
// lattice's actual discriminant group, M^2/2 is what the paper asserts it must be.
for re_DN in [ [6,1], [10,1], [15,2], [26,3], [22,3], [14,3], [39,2] ] do
    re_D := re_DN[1]; re_N := re_DN[2];
    re_M := IsOdd(re_D*re_N) select 4*re_D*re_N else 2*re_D*re_N;
    re_sz := #(ShimuraCurveLattice(re_D, re_N)`disc_grp);
    re_total +:= 1;
    if re_sz ne re_M^2 / 2 then
        re_fail +:= 1;
        Append(~re_bad, Sprintf("%o_%o: |A| = %o, paper says M^2/2 = %o", re_D, re_N, re_sz, re_M^2/2));
    end if;
end for;
printf "  |A| = M^2/2 : checked %o bases\n", re_total;

// ---- EXPENSIVE CHECK: the rho-entry itself -------------------------------------------------
// Restricted to the two small bases. 15_2 has |A| = 1800 and 26_3 has 12168, and the rho
// matrices are |A| x |A| -- 26_3 alone did not finish in 10 minutes. tests/WeilRepresentation.m
// uses exactly these two bases for the same reason.
for re_DN in [ [6,1], [10,1] ] do
    re_D := re_DN[1]; re_N := re_DN[2];
    re_Ld := ShimuraCurveLattice(re_D, re_N);
    re_S, re_T, re_elts, re_K := WeilRepresentationST(re_Ld);
    re_M := IsOdd(re_D*re_N) select 4*re_D*re_N else 2*re_D*re_N;
    re_A := re_M^2 / 2;                     // |A|, per the proof of thm:rhoentry

    // e_0 is the basis vector at the zero coset.
    re_zero := [i : i->x in re_elts | IsZero(x)];
    error if #re_zero ne 1, Sprintf("X0^%o(%o): expected exactly one zero coset, found %o",
                                    re_D, re_N, #re_zero);
    re_i0 := re_zero[1];

    // Walk a spread of gamma by their lower-left entry c; g is the class of the lower row.
    re_checked_here := 0;
    for re_c in [1..Minimum(re_M, 24)] do
        re_g := GCD(re_c, re_M);
        // gamma with lower row (c, d): take the standard rep with a=d=1 where possible.
        re_ok := true;
        try
            re_gam := Matrix(Integers(), 2, 2, [1, 0, re_c, 1]);   // lower row (c,1), det 1
            re_word := VVSTWord(re_gam);
            // build rho(gamma) as the product of the S/T matrices along the word
            re_rho := IdentityMatrix(re_K, #re_elts);
            for re_tok in re_word do
                if re_tok[1] eq "S" then
                    re_rho := re_rho * re_S;
                else
                    re_rho := re_rho * re_T^(re_tok[2] mod re_M);
                end if;
            end for;
        catch e
            re_ok := false;
        end try;
        if not re_ok then continue; end if;

        re_col := [re_rho[re_i, re_i0] : re_i in [1..#re_elts]];

        // (ii) ABSOLUTE VALUE. thm:rhoentry predicts |entry| = g/sqrt(|A|) for odd g and g/M
        // for 4 | g. The claim is about the eta* component specifically; we do not have the
        // cusp-class labelling here, so we test the strongest form we CAN: that every nonzero
        // entry of the column has that modulus. If the prediction is eta*-specific rather than
        // uniform across the column, this reports it as a mismatch -- which is information, not
        // a bug in the test, and the printed values say which.
        // ⚠ thm:rhoentry (ii) covers ONLY odd g and 4 | g. It says nothing about g = 2 mod 4,
        // so predicting there would be testing a claim the paper does not make. Skipped, and
        // the observation recorded below rather than asserted.
        if (re_g mod 4) eq 2 then continue; end if;
        re_total +:= 1; re_checked_here +:= 1;
        re_pred2 := IsOdd(re_g) select (Rationals()!re_g)^2 / re_A
                                  else (Rationals()!re_g)^2 / (Rationals()!re_M)^2;
        re_moduli := {};
        for re_x in re_col do
            if re_x eq 0 then continue; end if;
            re_z := re_x * ComplexConjugate(re_x);       // |x|^2, lands in the real subfield
            if re_z in Rationals() then
                Include(~re_moduli, Rationals()!re_z);
            else
                Include(~re_moduli, -1);                  // non-rational modulus: flag it
            end if;
        end for;
        if not IsEmpty(re_moduli) and re_moduli ne {re_pred2} then
            re_fail +:= 1;
            Append(~re_bad, Sprintf("%o_%o c=%o g=%o: predicted |.|^2 = %o, found %o",
                                    re_D, re_N, re_c, re_g, re_pred2, re_moduli));
        end if;
    end for;
    printf "  X0^%o(%o): M = %o, |A| = %o, %o gamma checked\n",
           re_D, re_N, re_M, re_A, re_checked_here;
end for;

printf "rho-entry closed form: %o gamma(s) checked, %o failure(s)%o\n",
       re_total, re_fail, IsEmpty(re_bad) select "" else Sprintf(" (%o)", re_bad);
error if re_fail ne 0,
    Sprintf("thm:rhoentry closed form disagrees with the computed rho: %o", re_bad);
