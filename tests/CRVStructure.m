// tests/CRVStructure.m
//
// THE FIRST CHECK OF ANY KIND ON `CRV` ENTRIES.
//
// ⚠ `VerifyModelSet` skips them outright --
//     if Type(e[2]) eq MonStgElt then continue; end if;   // "CRV": non-hyperelliptic entry
// -- so for every paired presentation in `data/models/` the stored equations have never been
// validated against anything, not even for internal consistency. 21 entries across 16 files.
//
// WHAT THIS CHECKS, and why it is cheap enough to be checked. Two structural facts that need no
// Shimura input at all, only the stored strings:
//   [1] the equations must define an IRREDUCIBLE scheme -- a `CRV` entry is a curve;
//   [2] the ambient WEIGHTS are derivable (a fibre coordinate's weight is half the degree of its
//       own equation, every other variable has weight 1), and with them the curve's genus must
//       equal the genus recorded beside it.
// Two shapes occur: the PAIRED one (y^2 and x^2 against forms in s, z, in P(1,2,1,1)) and, since
// 2026-10-03, the FIBRE-PRODUCT one (y_1^2, ..., y_k^2 against forms in s, z; FibreProductCovers.m),
// for which the curve is built as the tower of double covers of the t-line.
// [2] doubles as the check that the weights ARE derivable, which is what `PROVENANCE.md` needs
// before `ModelRegen` and the `X0_*` helper can handle `CRV` entries at all: they are currently
// skipped only because "the model file does not record the ambient weights". It does not have to
// -- 16 of 21 entries reconstruct exactly. No data migration is needed, just this derivation.
//
// ⚠⚠ WHAT IT FOUND, AND THE BUG BEHIND IT (2026-09-06 -- found, root-caused and FIXED same day).
// Five entries stored the SAME equation twice, up to `y` <-> `x`:
//
//     10_3  [1,10]  y^2 + 7/20*s^2 - 43/20*s*z + 2*z^2   and   x^2 + (the identical form)
//
// If `y^2 = q` and `x^2 = q` then `y^2 = x^2`, so `(y-x)(y+x) = 0` and the scheme is REDUCIBLE
// (measured: `IsIrreducible` false for all five). They could not be the genus-1 curves recorded
// beside them.
//
// ROOT CAUSE, in `EquationsAbovePointlessConics` (`EquationsCovers.m`). It builds the fibre product
//     C := Curve(P3, [y^2 - eqn2, x^2 - eqn1]);
// taking `eqn2` from a cover whose hyperelliptic polynomial has degree `g+1` and `eqn1` from a
// conic. At `g = 1` the required degree `g+1 = 2` IS a conic's degree, so the conic itself passed
// the degree test and could be picked for BOTH roles -- each of the five stored its own parent
// conic twice (`10_3 [1,10]` against `[1,2,5,10]`, and so on).
// FIXED by requiring the two roles be filled by different covers (`c ne other_curve`); those
// covers now DEFER, which is the honest outcome, and the five stored entries were emptied to match
// what the fixed pipeline produces.
//
// ⚠ The failure was invisible for as long as it was because `VerifyModelSet` skips every `CRV`
// entry -- there was no check of any kind on these until this file.

printf "Checking CRV entries for irreducibility and derivable weights...";

// The five entries measured degenerate on 2026-09-06 (identical equations up to y <-> x).
// Keep this list and the header in sync; shrink it as they are repaired.
// ✅ EMPTIED 2026-09-06 once the root cause was fixed (EquationsCovers.m: the two roles
// must be filled by DIFFERENT covers). The list is now empty and must STAY empty --
// a new degenerate entry is a regression, not something to add here.
CRV_KNOWN_BAD := [ Strings() | ];

crv_files := Split(Pipe("ls data/models/*.m | xargs grep -l '\"CRV\"'", ""), "\n");
crv_n := 0; crv_bad := []; crv_expected := [];

for crv_f in crv_files do
    crv_base := Split(Split(crv_f, "/")[#Split(crv_f, "/")], ".")[1];
    crv_base := Substring(crv_base, 8, #crv_base - 7);          // "models_D_N" -> "D_N"
    P<x> := PolynomialRing(Rationals());
    crv_models := eval (Read(crv_f) cat "\nreturn models;");
    for crv_k in Keys(crv_models) do
        for crv_e in crv_models[crv_k] do
            if Type(crv_e[2]) ne MonStgElt then continue; end if;
            crv_g := crv_e[1]; crv_s := crv_e[3];
            crv_tag := Sprintf("%o:[%o]", crv_base,
                               &cat[Sprint(t) cat "," : t in Sort(SetToSequence({z : z in crv_k}))]);
            crv_tag := Substring(crv_tag, 1, #crv_tag-2) cat "]";
            crv_n +:= 1;

            crv_ok := true; crv_got := -1; crv_irr := false;
            if exists{st : st in crv_s | Regexp("y[0-9]", st)} then
                // FIBRE-PRODUCT shape (FibreProductCovers.m): equations y_i^2 - F_i(s,z), one per
                // factor, in coordinates s, z of weight 1 and y_i of weight half the degree of its
                // own equation.  The curve is the compositum of the double covers y_i^2 = F_i(t,1)
                // of the t-line; building that tower IS the irreducibility check (Magma refuses a
                // reducible extension), and its genus must be the recorded one.
                R6<s, z, y1, y2, y3, y4> := PolynomialRing(Rationals(), 6);
                Pt<t> := PolynomialRing(Rationals());
                try
                    FF := RationalFunctionField(Rationals()); K := FF;
                    for i in [1..#crv_s] do
                        q := eval ("return " cat crv_s[i] cat ";");
                        F := R6.(2+i)^2 - q;                       // = F_i(s,z)
                        assert Evaluate(F, [R6.1, R6.2, 0, 0, 0, 0]) eq F;   // involves s, z only
                        fi := Evaluate(F, [t, 1, 0, 0, 0, 0]);
                        RK<Y> := PolynomialRing(K);
                        K := FunctionField(Y^2 - K!Evaluate(fi, FF.1));
                    end for;
                    crv_irr := true;
                    crv_got := Genus(K);
                catch e crv_ok := false;
                end try;
            else
                // PAIRED shape: y^2 and x^2 against forms in s, z; derive the weights from the
                // y-equation's degree
                R<yy,xx,ss,zz> := PolynomialRing(Rationals(), 4);
                crv_dy := Degree(eval ("return " cat crv_s[1] cat ";"))
                          where y is yy where x is xx where s is ss where z is zz;
                crv_wy := crv_dy div 2;
                try
                    Pw<xw,yw,sw,zw> := WeightedProjectiveSpace(Rationals(), [1,crv_wy,1,1]);
                    crv_eqs := [eval ("return " cat st cat ";") : st in crv_s]
                               where y is yw where x is xw where s is sw where z is zw;
                    crv_irr := IsIrreducible(Scheme(Pw, crv_eqs));
                    if crv_irr then crv_got := Genus(Curve(Pw, crv_eqs)); end if;
                catch e crv_ok := false;
                end try;
            end if;

            if crv_ok and crv_irr and (crv_got eq crv_g) then continue; end if;
            crv_why := (not crv_irr) select "scheme is REDUCIBLE"
                       else (crv_ok select Sprintf("genus %o, recorded %o", crv_got, crv_g)
                                    else "could not be built");
            if crv_tag in CRV_KNOWN_BAD then
                Append(~crv_expected, crv_tag cat " (" cat crv_why cat ")");
            else
                Append(~crv_bad, crv_tag cat ": " cat crv_why);
            end if;
        end for;
    end for;
end for;

// A test that checked nothing would also print a pass; say how many entries were examined.
error if crv_n eq 0,
    "CRVStructure: NO EVIDENCE -- found zero CRV entries to check, so nothing was verified.";
error if not IsEmpty(crv_bad),
    Sprintf("CRVStructure: %o CRV entry/entries are not irreducible curves of their recorded "
            * "genus, and are NOT in CRV_KNOWN_BAD: %o", #crv_bad, crv_bad);
error if #crv_expected ne #CRV_KNOWN_BAD,
    Sprintf("CRVStructure: CRV_KNOWN_BAD lists %o entries but %o were found degenerate (%o). If one "
            * "was REPAIRED, remove it from the list; the list must not go stale.",
            #CRV_KNOWN_BAD, #crv_expected, crv_expected);

printf " ok (%o CRV entries; %o reconstruct from derived weights, %o known-degenerate)\n",
       crv_n, crv_n - #crv_expected, #crv_expected;
