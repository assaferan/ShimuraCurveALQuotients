// Obstruction screen for EITHER parity: drives BorcherdsForms with DeficitScreen, which reports
// the deficit from INSIDE the intrinsic and so reaches the odd-D 0-side rows that deficit.m omits.
//
//   magma -b DD:=15 NN:=1 vvdata/weyl-campaign/deficit_odd.m < /dev/null
//
// ⚠⚠ READ `wdef`, NOT `deficit`.  `deficit` = Ncols - Rank asks whether EVERY vector is in the
// image; the search only asks it of a target supported on the CM points' coordinates.  Full column
// rank is sufficient, NOT necessary, and on odd D the gap decides the verdict: 55_1 reads
// deficit 3 / wdef 0 and BUILDS; 21_2 sits at deficit 2 across its whole ladder and BUILDS.
//
// ⚠ ODD D IS A DIFFERENT LADDER -- over m, not over the pole order P.  At fixed m the deficit
// GROWS with P (15_1, m = -3: 1 3 6 10 14 19 across P = 10..266), because the 0-side row block is
// fixed by m_choice while a deeper P keeps adding columns.  That growth is what the old
// "odd-D overestimate ~20" was really measuring.
//
// ⚠ On odd D this is fast only when it CLEARS: an obstructed verdict needs the whole m ladder, and
// the 0-side basis at deep m costs (pole order -D0*m, e.g. 4005 at 15_1).  And no odd-D base is
// confirmed obstructed by a real run, so an odd "obstructed" line means NOT CLEARED, nothing more.
//
// Known values it reproduces (check these before trusting a new number):
//   38_5 obstructed deficit 1 every rung | 34_3 clear 0 | 146_1 1 then 0
//   142_1, 158_1 wdef = deficit = 1 (both CONFIRMED obstructed by a real run)
//   15_1 clear at m = -7 (16 s), 55_1 clear at m = -15 -- odd bases that build
SetQuitOnError(true);
SetColumns(0);
AttachSpec("ShimuraQuotients.spec");
D := 15; N := 1;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;
Prec := 100;
if assigned PP then Prec := StringToInteger(PP); end if;
printf "DEFICIT BASE %o %o\n", D, N;
t0 := Realtime();
curves := GetHyperellipticCandidates();
Xstar := rep{X : X in curves | X`D eq D and X`N eq N and IsStarCurve(X)};
_ := BorcherdsForms(Xstar, curves : Prec := Prec, DeficitScreen);
printf "DEFICIT TOTAL %os\n", Realtime() - t0;
printf "DONE\n";
quit;
