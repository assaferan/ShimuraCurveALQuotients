# The genus-four quotient X_0^6(403)/<w_2, w_13, w_93>: a model derived from the star

`Y = X_0^6(403)/<w_2, w_13, w_93>` has genus four and is a double cover of the genus-two star
`C = X_0^6(403)*`, branched at the images of the fixed points of `w_3` and `w_403` (the CM points of
discriminants −3 and −403, one rational point each; the −3 points are the elliptic points of order
three). `tests/Hyperelliptic6_403.m` proves `Y` hyperelliptic, so the branch divisor is a fibre of
`C`'s hyperelliptic map; `C` has one rational Weierstrass point (`x = 1`) and no rational torsion,
so a single normal form `y^2 = q1(t^2/c)`, `c` a twist supported on the bad primes. The star's
equation is the one fitted to its differentials (`star-model.log`, 58 good primes against the
newform), a candidate; the genus-four model inherits that status. Run from the repository root:

| script | what it shows |
|---|---|
| `enumerate.m` | the generic enumeration (every rational Weierstrass base, every fibre of height ≤ 10^6, every twist supported on {−1, 2, 3, 13, 31}): `magma -b DD:=6 NN:=403 FC:="[9,-21,33,-47,15,9,1]" HC:="[0,1,1]" HH:="{1,2,13,26,93,186,1209,2418}" tools/models-6_403/enumerate.m` — exactly one triple has the counts of `Y`: base `x = 1`, fibre `x = −7`, `c = 1` |
| `thuemahler.m` | the branch fibre, exactly: the weighted y-coordinate of the fibre is a unit outside {2, 3, 13, 31} and the quintic factor of the star's sextic is an S-unit (the resultant with `x − 1` is 72), a Thue–Mahler equation whose 17 solutions (`thuemahler.sols`, Gherga–Siksek's solver, 293 s) give 6 possible fibres, of which only `x = −7` passes the count test |
| `verify.m` | the resulting curve in its nice form, its two genus-two quotients (the star and `C'`), its counts against the trace formula at the 31 good primes below 150, and the conductor of `C'` |
