// FEASIBILITY PROBE for extending deficit.m to ODD D.
//
//   magma -b DD:=15 NN:=1 vvdata/weyl-campaign/oddprobe.m < /dev/null
//
// PLAN.md says the odd-D extension is "extract BorcherdsForms.m:876-978 and import it".  Reading
// the block shows it needs THREE things the even-D screen never computes:
//   * a SECOND WeaklyHolomorphicBasis call, the Zero one (nE0, eta_quotients_oo);
//   * `pts` from RationalandQuadraticCMPoints -- the field-of-definition work the screen exists
//     to skip -- because all_ms, hence m_choice, hence the 0-side pole_order, comes from it;
//   * therefore a sweep over (P, m_choice) PAIRS, not over P alone.
// This script measures what those cost before anything is refactored.
SetQuitOnError(true);
SetColumns(0);
AttachSpec("ShimuraQuotients.spec");
_ := ClassNumberLU(-4);          // force ClassNumberData.m to load before the import (CLAUDE.md)
import "BorcherdsForms.m" : get_D0_M_g;
D := 15; N := 1;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;
Prec := 100;
printf "ODDPROBE BASE %o %o  (M = %o)\n", D, N, 2*D*N;

curves := GetHyperellipticCandidates();
Xstar := rep{X : X in curves | X`D eq D and X`N eq N and IsStarCurve(X)};

t := Realtime();
E, n, n0, tt, eta_quotients := WeaklyHolomorphicBasis(D, N : Prec := Prec);
printf "ODDPROBE WHB_oo %os  (n = %o, n0 = %o, #eta = %o)\n", Realtime()-t, n, n0, #eta_quotients;

t := Realtime();
E0, nE0, _, eta_quotients_oo, eta_quotients_0 := WeaklyHolomorphicBasis(D, N : Prec := Prec, Zero, n0 := n0);
printf "ODDPROBE WHB_0  %os  (nE0 = %o, #eta_oo = %o, #eta_0 = %o)\n",
       Realtime()-t, nE0, #eta_quotients_oo, #eta_quotients_0;

t := Realtime();
pts, _ := RationalandQuadraticCMPoints(Xstar : bd := 2);
printf "ODDPROBE CMPOINTS %os  (#pts = %o)\n", Realtime()-t, #pts;
if #pts lt 3 then
    t := Realtime();
    pts, _ := RationalandQuadraticCMPoints(Xstar : bd := 2, coprime_to_level := false);
    printf "ODDPROBE CMPOINTS_noncoprime %os  (#pts = %o)\n", Realtime()-t, #pts;
end if;

all_ms := &cat[[(d[1] mod 4 eq 0) select d[1] div 4 else d[1] : d in pts] : pt in pts];
all_ms := Reverse(Sort([m : m in Set(all_ms)]));
printf "ODDPROBE ALL_MS %o\n", all_ms;
D0 := get_D0_M_g(D, N);
printf "ODDPROBE D0 %o  zero-side pole orders -D0*m = %o\n", D0, [-D0*m : m in all_ms];
k := -Valuation(qExpansionAtoo(tt, 1));
printf "ODDPROBE k %o  floor_pole %o\n", k, n0 + k - 1;
printf "DONE\n";
quit;
