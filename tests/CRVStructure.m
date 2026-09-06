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
//   [1] the two equations must define an IRREDUCIBLE scheme -- a `CRV` entry is a curve;
//   [2] the ambient WEIGHTS are derivable (`y`'s weight is half the degree of its own equation,
//       every other variable has weight 1), and with them the curve's genus must equal the genus
//       recorded beside it.
// [2] doubles as the check that the weights ARE derivable, which is what `PROVENANCE.md` needs
// before `ModelRegen` and the `X0_*` helper can handle `CRV` entries at all: they are currently
// skipped only because "the model file does not record the ambient weights". It does not have to
// -- 16 of 21 entries reconstruct exactly. No data migration is needed, just this derivation.
//
// ⚠⚠ WHAT IT FOUND, 2026-09-06. FIVE entries store the SAME equation twice, up to `y` <-> `x`:
//
//     10_3  [1,10]  y^2 + 7/20*s^2 - 43/20*s*z + 2*z^2   and   x^2 + (the identical form)
//
// If `y^2 = q` and `x^2 = q` then `y^2 = x^2`, so `(y-x)(y+x) = 0` and the scheme is REDUCIBLE --
// measured, `IsIrreducible` is false for all five. They cannot be the genus-1 curves recorded
// beside them. A fibre product of a double cover with ITSELF is reducible by construction, so this
// looks like two covers with the same equation being paired as if independent.
// These are listed in CRV_KNOWN_BAD so the test documents them rather than failing CI, exactly as
// `ModelRegen`'s MR_KNOWN_DRIFT does. ⇒ They are a REAL DEFECT to repair, not a convention to
// accommodate; when one is fixed, remove it from the list.

printf "Checking CRV entries for irreducibility and derivable weights...";

// The five entries measured degenerate on 2026-09-06 (identical equations up to y <-> x).
// Keep this list and the header in sync; shrink it as they are repaired.
CRV_KNOWN_BAD := [ "10_3:[1,10]", "10_3:[1,15]", "10_3:[1]", "10_3:[1,6]", "22_3:[1,3]" ];

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

            // derive the weights from the y-equation's degree
            R<yy,xx,ss,zz> := PolynomialRing(Rationals(), 4);
            crv_dy := Degree(eval ("return " cat crv_s[1] cat ";"))
                      where y is yy where x is xx where s is ss where z is zz;
            crv_wy := crv_dy div 2;

            crv_ok := true; crv_got := -1; crv_irr := false;
            try
                Pw<xw,yw,sw,zw> := WeightedProjectiveSpace(Rationals(), [1,crv_wy,1,1]);
                crv_eqs := [eval ("return " cat st cat ";") : st in crv_s]
                           where y is yw where x is xw where s is sw where z is zw;
                crv_irr := IsIrreducible(Scheme(Pw, crv_eqs));
                if crv_irr then crv_got := Genus(Curve(Pw, crv_eqs)); end if;
            catch e crv_ok := false;
            end try;

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
