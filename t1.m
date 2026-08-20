// Step 1 check: coset transversal + ST words.
AttachSpec("ShimuraQuotients.spec");
Attach("vvlib.m");

for M in [12, 20, 60, 84] do
    reps := VVCosetReps(M);
    idx := Index(Gamma0(M));
    // cosets must be pairwise distinct: g1 g2^-1 in Gamma_0(M) iff M | lower-left
    bad := 0;
    for i in [1..#reps] do
        for j in [i+1..#reps] do
            h := reps[i]*reps[j]^(-1);
            if h[2][1] mod M eq 0 then bad +:= 1; end if;
        end for;
    end for;
    // words reconstruct
    wbad := 0; slens := [];
    for g in reps do
        w := VVSTWord(g);
        if VVWordMatrix(w) ne g then wbad +:= 1; end if;
        Append(~slens, #[t : t in w | t[1] eq "S"]);
    end for;
    printf "M=%-4o  #reps=%-5o Index(Gamma0)=%-5o  collisions=%-3o  word_fail=%-3o  maxS=%o\n",
           M, #reps, idx, bad, wbad, Maximum(slens);
end for;
quit;
