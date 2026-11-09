// tests/CRVTraceCounts.m
//
// EVERY curve stored as a `CRV` entry, against the Eichler-Selberg trace formula.
//
// The non-hyperelliptic entries of `data/models/` -- the older pairs over a conic and the fibre
// products of 2026-10-03/04 -- were accepted, when they were built, on the genus and on point
// counts at a few primes, and nothing afterwards looked at them again.  `VerifyModelSet` skips
// them, and `tests/CRVStructure.m` only asks that the stored equations define an irreducible curve
// of the recorded genus.  A factor rewritten in another coordinate on the star curve, or an
// equation edited by hand, gives a curve of the right genus that is not the right curve, and until
// now nothing in the suite would have noticed: the committed strings had never been compared with
// anything outside the construction that produced them.
//
// WHAT IS CHECKED.  A `CRV` entry is the fibre product over the star line of the double covers
// y_i^2 = f_i(t) its equations record (`CRVDoubleCovers` reads them back from the stored strings
// alone).  For each entry this test builds that tower over F_p at the first two good primes and
// compares the number of places of degree 1, and of degree 1 and 2 together, with the number of
// points of X_0(D,N)/W over F_p and F_{p^2} that the Eichler-Selberg trace formula gives on the
// W-fixed part of the D-new space (`ComputePointsViaTrace`).  That formula is evaluated from D, N
// and W and knows nothing about any model, so it is the independent side of the comparison; the
// genus is the Shimura-curve genus formula's.
//
// A prime dividing DN is not used, and a prime at which the tower cannot be built or the genus
// drops is bad reduction of this model, not a verdict: it is replaced by the next prime, searching
// the first 40.  Measured on 2026-10-13: 299 entries, 598 (entry, prime) pairs, 147 s on a Mac, so
// every entry gets its two primes.
//
// ⚠ TWO PRIMES IS A SCREEN, NOT A PROOF, and the screen does let things through.  Twisting one
// factor of the stored X_0^6(23)/w_23 by -1 is caught at p = 7, but scaling it by 2 agrees at
// p = 5 and p = 7 and is only caught at p = 11 (measured 2026-10-13).  What the test is for is a
// factor rewritten in another coordinate or an equation edited by hand, which these two primes do
// reject; a planted quadratic twist needs more of them.

printf "Checking CRV entries against the trace formula...";

ctc_curves := GetHyperellipticCandidates();
ctc_files := Split(Pipe("ls data/models/*.m | xargs grep -l '\"CRV\"'", ""), "\n");
ctc_nent := 0; ctc_npair := 0; ctc_bad := []; ctc_thin := [];

for ctc_f in ctc_files do
    if #ctc_f eq 0 then continue; end if;
    ctc_parts := Split(Split(ctc_f, "/")[#Split(ctc_f, "/")], "_");
    ctc_D := StringToInteger(ctc_parts[2]);
    ctc_N := StringToInteger(Split(ctc_parts[3], ".")[1]);
    P<x> := PolynomialRing(Rationals());
    ctc_models := eval (Read(ctc_f) cat "\nreturn models;");
    for ctc_k in Keys(ctc_models) do
        for ctc_e in ctc_models[ctc_k] do
            if Type(ctc_e[2]) ne MonStgElt then continue; end if;
            ctc_g := ctc_e[1];
            ctc_W := {Integers()| a : a in ctc_k};
            ctc_tag := Sprintf("%o_%o W=%o", ctc_D, ctc_N, Sort(SetToSequence(ctc_W)));
            ctc_nent +:= 1;
            ctc_fs := CRVDoubleCovers(ctc_e[3]);
            ctc_X := [Y : Y in ctc_curves | Y`D eq ctc_D and Y`N eq ctc_N and Y`W eq ctc_W];
            error if #ctc_X ne 1,
                Sprintf("CRVTraceCounts: %o matches %o curve records", ctc_tag, #ctc_X);
            ctc_X := ctc_X[1];
            error if ctc_X`g ne ctc_g,
                Sprintf("CRVTraceCounts: %o is stored with genus %o, the genus formula gives %o",
                        ctc_tag, ctc_g, ctc_X`g);
            ctc_done := 0; ctc_p := 2; ctc_tried := 0;
            while ctc_done lt 2 and ctc_tried lt 40 do
                ctc_p := NextPrime(ctc_p); ctc_tried +:= 1;
                if (ctc_D * ctc_N) mod ctc_p eq 0 then continue; end if;
                ctc_Kp := RationalFunctionField(GF(ctc_p)); ctc_L := ctc_Kp;
                ctc_built := true;
                try
                    for ctc_fi in ctc_fs do
                        RY<Y> := PolynomialRing(ctc_L);
                        ctc_L := FunctionField(Y^2 - ctc_L!Evaluate(
                                     PolynomialRing(GF(ctc_p))!ctc_fi, ctc_Kp.1));
                    end for;
                catch err ctc_built := false;
                end try;
                if not ctc_built then continue; end if;
                if Genus(ctc_L) ne ctc_g then continue; end if;
                ctc_cnt := [&+[e * #Places(ctc_L, e) : e in Divisors(d)] : d in [1..2]];
                ctc_exp := [ComputePointsViaTrace(ctc_X, ctc_p, d) : d in [1..2]];
                if ctc_cnt ne ctc_exp then
                    Append(~ctc_bad, Sprintf("%o, p = %o: the curve has %o points over F_p, F_p^2, "
                                             * "the trace formula %o", ctc_tag, ctc_p, ctc_cnt, ctc_exp));
                end if;
                ctc_done +:= 1; ctc_npair +:= 1;
            end while;
            if ctc_done lt 2 then
                Append(~ctc_thin, Sprintf("%o (%o good primes)", ctc_tag, ctc_done));
            end if;
        end for;
    end for;
end for;

// A test that checked nothing would also print a pass; say how much was examined.
error if ctc_nent eq 0 or ctc_npair eq 0,
    Sprintf("CRVTraceCounts: NO EVIDENCE -- %o CRV entries and %o point-count comparisons, so "
            * "nothing was verified.", ctc_nent, ctc_npair);
error if not IsEmpty(ctc_bad),
    Sprintf("CRVTraceCounts: %o stored CRV curve(s) disagree with the trace formula: %o",
            #ctc_bad, ctc_bad);
error if not IsEmpty(ctc_thin),
    Sprintf("CRVTraceCounts: %o entry/entries had fewer than two good primes among the first 40 "
            * "(so they were barely checked): %o", #ctc_thin, ctc_thin);

printf " ok (%o entries, %o (entry, prime) pairs, F_p and F_p^2 each)\n", ctc_nent, ctc_npair;
