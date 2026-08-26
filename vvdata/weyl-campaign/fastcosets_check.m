// Validation gate for fastcosets.m: fastCosetReps must equal VVCosetReps as a
// SEQUENCE OF MATRICES (not merely as a set of cosets) at every level the
// campaign has used, so it is a drop-in with no gauge change at all.
//   magma -b vvdata/weyl-campaign/fastcosets_check.m
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/fastcosets.m";

allok := true;
for M in [60, 84, 132, 140, 420, 660] do
    t := Cputime(); slow := VVCosetReps(M); ts := Cputime(t);
    t := Cputime(); fast := fastCosetReps(M); tf := Cputime(t);
    same := (#slow eq #fast) and forall{ i : i in [1..#slow] | slow[i] eq fast[i] };
    if not same then allok := false; end if;
    printf "M = %4o: %5o cosets  VVCosetReps %8o s   fastCosetReps %8o s   %o\n",
        M, #fast, ts, tf, same select "IDENTICAL" else "*** DIFFER ***";
end for;
printf "FASTCOSETS VERDICT: %o\n", allok select "PASS" else "FAIL";

// and the target level, which the slow routine cannot reach
t := Cputime(); f := fastCosetReps(4620); tf := Cputime(t);
printf "M = 4620: %o cosets in %o s (expected 13824)\n", #f, tf;
error if #f ne 13824, "wrong coset count at 4620";
// spot-check that they really are distinct cosets of Gamma_0(4620)
M := 4620;
keys := { };
for g in f do
    c := g[2][1] mod M; d := g[2][2] mod M;
    orb := { <(u*c) mod M, (u*d) mod M> : u in [1..M] | GCD(u, M) eq 1 };
    Include(~keys, Min([ x : x in orb ]));
end for;
printf "M = 4620: %o distinct P^1(Z/M) classes among %o reps  %o\n",
    #keys, #f, (#keys eq #f) select "OK" else "*** COLLISION ***";
error if #keys ne #f, "reps are not pairwise inequivalent";
printf "DONE\n";
quit;
