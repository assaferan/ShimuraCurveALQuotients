// tests/_crossnorm.m -- NOT a test (leading underscore: excluded from the CI matrix). Run by hand:
//     magma -b Dd:=10 Nn:=3 tests/_crossnorm.m < /dev/null
//
// ⇒ ARE THE CM VALUES RIGHT, WHEN THERE IS NO PUBLISHED VALUE TO COMPARE AGAINST?
//
// This answers that WITHOUT an external oracle, by running the pipeline in two trees that differ
// only in which (infty, P, Q) anchor the Borcherds sweep picks, dumping the hauptmodul's values at
// every CM point in the Schofer table, and checking that ONE Mobius map carries one set onto the
// other. Written 2026-09-24 to settle whether 0ca6e37's newly admitted CM points carry WRONG
// VALUES or merely a DIFFERENT NORMALISATION. Answer at 10_3: the same function, so the values are
// right and the cover loss is a normalisation artefact. See HANDOFF 2026-09-24.
//
// ⚠⚠ THE OBVIOUS FORM OF THIS TEST IS VACUOUS, AND THE FIRST DRAFT OF IT WAS.
// Any two hauptmoduls of a genus-0 curve are Mobius-related BY DEFINITION, so "does a map exist"
// proves nothing, and fitting on the three shared anchors determines the map exactly and checks
// nothing. It bites only because the Schofer table is evaluated at MORE discriminants than the
// three the map needs: fit on three, and every remaining one is a genuine, over-determined check.
// ⇒ If #Discs <= 3 this tool CANNOT conclude anything, and it says so rather than printing a
// reassuring nothing.
//
// USAGE, which needs two trees because the anchor choice lives in the SOURCE (BorcherdsForms.m,
// the `pts` ordering), not in a flag:
//     1. run it in tree A, keep the VAL lines;
//     2. run it in tree B, keep the VAL lines;
//     3. fit the Mobius map on any three shared discriminants and check the rest.
// The fitting step is deliberately left outside: the two runs are minutes apart in different trees,
// and a tool that silently fitted whatever it found would hide a disc-set mismatch. Compare the
// disc sets FIRST -- if they differ, the comparison is between different objects.
//
// ⚠ Cross-ratio convention, matching tests/ExternalCMValues.m: the map sending z0 -> 0, z1 -> oo,
// z2 -> 1 is  (z - z0)/(z - z1) * (z2 - z1)/(z2 - z0),  with each factor dropped when its argument
// is Infinity.
AttachSpec("ShimuraQuotients.spec");
D := StringToInteger(Dd); N := StringToInteger(Nn);
curves := GetHyperellipticCandidates();
assert exists(Xstar){X : X in curves | X`D eq D and X`N eq N and IsStarCurve(X)};

tab  := ValuesAtCMPoints(Xstar, curves : Prec := 100);
ds   := tab`Discs;
si   := tab`sIndex;
vals := tab`Values;

printf "X0^%o(%o)  hauptmodul row %o,  %o discriminant(s) in the table\n", D, N, si, #ds;
if #ds le 3 then
    printf "⚠ ONLY %o DISCRIMINANT(S): three fit the Mobius map exactly, so there is nothing left to\n", #ds;
    printf "  check and this run CANNOT support any conclusion. Raise MaxNum or pick another base.\n";
else
    printf "⇒ %o over-determined check(s) available once three are spent fitting the map.\n", #ds - 3;
end if;
for j -> d in ds do
    printf "VAL %o %o\n", d, vals[si][j];
end for;
exit;
