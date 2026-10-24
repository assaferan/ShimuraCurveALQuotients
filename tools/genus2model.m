// A MODEL OVER Q FOR A GENUS-TWO ATKIN-LEHNER QUOTIENT, FROM ITS DIFFERENTIALS.
//
//   magma -b DD:=2002 NN:=1 WS:=1,2,7,11,13,14,22,26,77,91,143,154,182,286,1001,2002 \
//         CUTA:=60 PRA:=60 tools/genus2model.m < /dev/null
//
// A thin driver over NelsonModelFromDifferentials, which does the work and is where the method is
// documented; this prints its diagnostics in a readable form.  Optional: SAMPLES (default 10),
// RATDISCS to override the CM discriminants the intrinsic would find for itself (a comma-separated
// list; normalising on a point that is NOT rational yields a plausible wrong curve that nothing
// downstream detects, which is why the default is to let it look), PRIMECAP for the twist (300).
// The method is written up in paper/nelson-method.tex, "A genus-two quotient, and a model over Q".

AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
D := StringToInteger(DD); N := StringToInteger(NN);
W := { StringToInteger(x) : x in Split(WS, ",") };
CUT := StringToInteger(CUTA); PR := StringToInteger(PRA);
ns := assigned SAMPLES select StringToInteger(SAMPLES) else 0;
rds := assigned RATDISCS select [ StringToInteger(x) : x in Split(RATDISCS, ",") ] else [];
pc := assigned PRIMECAP select StringToInteger(PRIMECAP) else 300;
R8 := RealField(8); C12 := ComplexField(12);

printf "=== X_0^%o(%o)/%o, cutoff %o, %o digits\n", D, N, Sort(Setseq(W)), CUT, PR;
C, dg := NelsonModelFromDifferentials(D, N, W : Precision := PR, Cutoff := CUT, Samples := ns,
                                               Discriminants := rds, PrimeCap := pc);
printf "\nrational CM points used, with their coordinates:\n";
for i in [1..#dg`discs] do
    printf "   d = %-7o %-10o x = %o\n", dg`discs[i],
           dg`isbranch[i] select "fixed pt" else "ordinary", C12 ! dg`cmx[i];
end for;
printf "spread of Im(x)/Re(x) over them: %o   (a rational point must lie on that one line)\n",
       R8 ! dg`spread;
printf "worst prediction at the points held back from the fit: %o\n", R8 ! dg`worst;
printf "tolerance %o, denominator bound %o; margins over a chance fit: %o\n",
       R8 ! dg`tol, dg`bound, [ R8 ! m : m in dg`margins ];
printf "\nsextic in the coordinate fixed by the CM points: %o\n", dg`sextic;
printf "absolute invariants: %o\n", dg`invariants;
printf "twisting constant %o; counts against the newform: %o good primes, %o mismatches\n",
       dg`twist, dg`counts[1], dg`counts[2];
printf "\nTHE MODEL: %o\n", C;
quit;
