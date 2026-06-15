// Direct measurement of the class-number-table speedup on the primitive the pipeline uses
// (ClassNumberBatchLU), across |d| scales.  Run twice (env CLASS_GROUPS_DIR off/on).
// SetSeed fixes the discriminant sample so both runs request identical data.
SetQuitOnError(true);
AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 0);
dir := GetEnv("CLASS_GROUPS_DIR");
printf "=== MODE: %o ===\n", (dir ne "") select dir else "no tables (live)";

NS := 20000;
for scale in [10^5, 10^7, 10^8, 2^28] do
    SetSeed(12345);                       // identical sample per scale across both runs
    Ds := [];
    while #Ds lt NS do
        d := Random(100, scale);
        if (d mod 4) in [0,3] then Append(~Ds, -d); end if;   // -d = 0,1 mod 4
    end while;
    t := Realtime();
    res := ClassNumberBatchLU(Ds);
    dt := Realtime() - t;
    printf "|d| <= %-12o : %o discriminants, %8.3o s  (%.4o ms/disc)\n",
        scale, #Ds, dt, 1000*dt/#Ds;
end for;
quit;
