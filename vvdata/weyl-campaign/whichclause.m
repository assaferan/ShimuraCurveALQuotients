// gtsweep's "NOT MEASURABLE" has TWO causes and the message does not say which fired:
//   (a) fewer than two INFORMATIVE points  (SchoferFormula.m:1272, early return {})
//   (b) no (r1,r2) in [-12..12]^2 fits the table -- a genuine INCONSISTENCY.
// (a) is a limitation of the sweep; (b) is a defect of the base.  Replay the logged tables
// through the same code path and report which one it is.  No BorcherdsForms run needed.
AttachSpec("ShimuraQuotients.spec");

// D, N, ds, s, stilde, degs -- transcribed from the gtsweep logs.
cases := [* *];
Append(~cases, <33, 2, [-4,-12,-88,-132,-148,-168,-15,-55],
    [* Infinity(), 0, 16/9, 3, 1, 1, 4, 16 *],
    [* Infinity(), 64, 0, 11, 25, 7, 20, 176 *], [1,1,1,1,1,1,1,2]>);
Append(~cases, <46, 3, [-3,-8,-24,-123,-312,-372,-35,-35],
    [* 0, 3/8, 1/2, Infinity(), 3/7, 1/4, 9/28, 9/28 *],
    [* 1, 0, 1/3, Infinity(), 1/7, 1/3, 2/7, 2/7 *], [1,1,1,1,1,1,2,2]>);
Append(~cases, <35, 2, [-8,-280,-88,-120,-148,-168,-7,-15],
    [* Infinity(), 0, 11/5, 1, 5, 3/5, 7/5, 3 *],
    [* Infinity(), 11/4, 0, 4, 7/2, 2, 1, 1 *], [1,1,1,1,1,1,1,1]>);
// 33_2 AFTER the integrality fix, run against the MAIN package (which has
// M0MultiplierExact -- the campaign branch does not).  Before the fix this base was
// clause (b), genuinely inconsistent; the question is whether the fix moved it to (a).
Append(~cases, <33, 2, [-4,-88,-132,-148,-168,-232,-15,-55],
    [* Infinity(), 11/9, 0, 4, 2, 2/9, 1, 11 *],
    [* Infinity(), 25/9, 4, 0, 2, 34/9, 5, 29 *], [1,1,1,1,1,1,1,2]>);
Append(~cases, <21, 2, [-4,-84,-168,-88,-120,-148,-7,-15],
    [* Infinity(), 0, 9, 24, 8, 64/3, 12, 4 *],
    [* Infinity(), 9, 0, 33, 1, 37/3, 21, 5 *], [1,1,1,1,1,1,1,1]>);

for c in cases do
    D := c[1]; N := c[2]; ds := c[3]; s := c[4]; stilde := c[5]; degs := c[6];
    i_s0 := Index(s, 0); i_st0 := Index(stilde, 0);
    scale := s[i_st0]; scale_tilde := stilde[i_s0];
    fires := func<d | not IsEmpty(PrimeDivisors(N div GCD(N, FundamentalDiscriminant(d))))>;
    kA := fires(ds[i_st0]) select 1 else 0;
    kB := fires(ds[i_s0])  select 1 else 0;
    rows := []; exps := []; which := [];
    for i->d in ds do
        if degs[i] ne 1 or i eq i_s0 or i eq i_st0 then continue; end if;
        if s[i] cmpeq Infinity() or stilde[i] cmpeq Infinity() then continue; end if;
        k := fires(d) select 1 else 0;
        Append(~rows, [Rationals() | s[i]/scale, stilde[i]/scale_tilde]);
        Append(~exps, [k - kA, k - kB]);
        Append(~which, d);
    end for;
    ninf := #[j : j in [1..#rows] | exps[j] ne [0,0]];
    printf "\n=== X0^%o(%o): %o usable rows, %o informative\n", D, N, #rows, ninf;
    // per-row: is THIS row satisfiable on its own, for some (r1,r2)?
    bad := [];
    for j := 1 to #rows do
        okj := exists{<r1,r2> : r1, r2 in [-12..12] |
                 exists{<g1,g2> : g1, g2 in [1,-1] |
                   g1*rows[j][1]*(Rationals()!N)^(r1*exps[j][1])
                 + g2*rows[j][2]*(Rationals()!N)^(r2*exps[j][2]) eq 1}};
        printf "   d = %-6o exps %o  row %o  satisfiable-alone: %o\n",
               which[j], exps[j], rows[j], okj;
        if not okj then Append(~bad, which[j]); end if;
    end for;
    sols := HauptmodulM0Residuals(s, stilde, ds, degs, N);
    cause := ninf lt 2 select "(a) TOO FEW INFORMATIVE POINTS -- a sweep limitation"
             else (IsEmpty(sols) select "(b) INCONSISTENT -- no multiplier fits: a DEFECT of the base"
                                 else "consistent");
    printf "   sols = %o  (#%o)\n   CAUSE %o %o : %o\n", sols, #sols, D, N, cause;
    if not IsEmpty(bad) then printf "   individually unsatisfiable rows: %o\n", bad; end if;
end for;
printf "\nDONE\n";
quit;
