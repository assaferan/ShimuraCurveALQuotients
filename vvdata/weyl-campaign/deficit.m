// PREDICTOR PROTOTYPE: compute the Borcherds-obstruction deficit WITHOUT the CM points and
// WITHOUT the divisor-triple loop.
//
//   magma -b DD:=38 NN:=5 deficit.m
//
// The deficit is  Ncols(mat) - Rank(ech_basis * mat)  and depends only on the weakly
// holomorphic basis and coeffs_to_divisor_matrix -- NOT on any divisor choice.  This script
// therefore skips RationalandQuadraticCMPoints (field-of-definition work) and the 96-triple
// search entirely.  It also times each stage, so we learn whether a cheap predictor is even
// possible: if WeaklyHolomorphicBasis dominates, the predictor costs what the pipeline costs.
//
// Ground truth from the spanprobe runs:
//   38_5 poleord 190 -> deficit 1     38_7 poleord 134/266 -> 0     34_3 poleord 51/102 -> 0
AttachSpec("ShimuraQuotients.spec");
// REPAIRED 2026-09-02.  This script used to call two helpers that existed ONLY in the
// -spanprobe worktree's patched BorcherdsForms.m.  That worktree was retired on 2026-08-30
// (content preserved as bfprof-instrumentation.patch), so the script had been dead since --
// and it died AFTER paying the full WeaklyHolomorphicBasis cost, 20+ minutes on a large base.
// The production functions are file-local, so import them directly rather than depending on a
// patch or a deleted worktree.
import "BorcherdsForms.m" : basis_of_weakly_holomorphic_forms, coeffs_to_divisor_matrix;
// VALIDITY: EVEN D ONLY.  For odd D, BorcherdsForms vertically joins a 0-side block
// (coeffs_0_oo) into coeffs_trunc; this script computes only the infinity part, so the extra
// rows are missing and the reported deficit is an OVERESTIMATE.  The symptom is that it drifts
// with the pole order instead of staying invariant -- 65_2 gives 5, 6, 6, 9 at P = 62, 102,
// 134, 190, where a genuine deficit is constant (38_5: 1 at every P).  Do not read odd-D
// numbers from this script until the 0-side block is included.
D := 38; N := 5;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;
Prec := 100;
if assigned PP then Prec := StringToInteger(PP); end if;

printf "DEFICIT BASE %o %o\n", D, N;

t0 := Realtime();
E, n, n0, t, eta_quotients := WeaklyHolomorphicBasis(D, N : Prec := Prec);
t_wh := Realtime() - t0;
printf "DEFICIT TIMING WeaklyHolomorphicBasis %os  (n = %o, n0 = %o, #eta = %o)\n",
       t_wh, n, n0, #eta_quotients;

k := -Valuation(qExpansionAtoo(t, 1));
printf "DEFICIT k %o\n", k;

// Sweep pole orders: the deficit should be stable, and we want it at the pole orders the real
// runs used.  min_m = -(pole order); the pipeline also floors min_m at -(n0 + k - 1).
floor_pole := n0 + k - 1;
poles := [floor_pole];
for P in [51, 102, 134, 190, 266] do
    if (P gt floor_pole) and (P notin poles) then Append(~poles, P); end if;
end for;
Sort(~poles);

for P in poles do
    tt := Realtime();
    try
        ech, etas, T := basis_of_weakly_holomorphic_forms(P, eta_quotients, n0 + 1, n, t);
        mat, rds := coeffs_to_divisor_matrix(-P, D, N, Ncols(ech));
        ct := ChangeRing(ech, Rationals()) * ChangeRing(mat, Rationals());
        r := Rank(ct);
        printf "DEFICIT P %o rows %o cols %o nds %o rank %o deficit %o  (%os)\n",
               P, Nrows(ct), Ncols(ct), #rds, r, Ncols(ct) - r, Realtime() - tt;
    catch e
        printf "DEFICIT P %o SKIP %o\n", P, e`Object;
    end try;
end for;
printf "DEFICIT TOTAL %os\n", Realtime() - t0;
printf "DONE\n";
quit;
