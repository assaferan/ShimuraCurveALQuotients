# The trace formula at gcd(n, N) > 1 -- probes behind PRs #56 and #57 (2026-09-30)

Run any script from a code tree: `magma -b <script>.m < /dev/null`.  Each .log is that script's
output on the tree named in its header (main @ 95cf87b unless a PR branch is named).

    oracle_scan.m     newform trace (formula) vs modular symbols, levels 12..108 -- the 72/336 grid;
                      oracle_before/after.log: identical on main and on the #56 branch
    popa_scan.m       splits the discrepancy: Popa full-space formula wrong on 80/208, newform 72/208
    cor55_oracle.m    Assaf Cor 5.5 run on MODULAR-SYMBOL full traces: 0/208 -> the corollary is right
    sfast_probe.m     Sfast vs brute-force S with p | (n,N): 50/154 differ; Popa with brute-force S: 0/208
    sfast_fix.m       the unit-root fix: 0/3440 vs brute force, 0/828 Popa vs modsym (levels 8..108, k 2,4,6)
    old_sfast_count.m the mismatch count of the OLD Sfast on tests/SfastUnitCondition.m's exact loop: 774/2677
    modsym_oracle.m   the 21 TraceDNewALFixed values from modular symbols (D-new, W-averaged, JL sign):
                      the 5 with gcd(n, DN) > 1 differed on main -- the D-new assembly kept only the
                      n' = 1 term of Cor 4.27
    dnew_general.m    the general-n D-new trace (Lemma 4.20 block by block, #57): 0/1548 differ at
                      gcd(n, DN) > 1 across 23 (D, N), several W, n in {2..25}; 0/242 coprime

⚠ Oracle trap: Magma's HeckeOperator on a SUBSPACE dies at (N,k,n) = (6,4,2); restrict the ambient
operator (Solution(B, B*T)) as tests/trace_formula.m does.
