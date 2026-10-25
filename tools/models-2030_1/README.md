# The genus-four quotient X_0^2030(1)/<w_10, w_14, w_58>: a model derived from the star

`Y = X_0^2030(1)/<w_10, w_14, w_58>` has genus four and is a double cover of the genus-two star
`C = X_0^2030(1)*`, branched at the images of the fixed points of `w_2` and `w_70` (the CM points of
discriminants −8 and −280, one rational point each). `tests/Hyperelliptic2030_1.m` proves `Y`
hyperelliptic, so the branch divisor is a fibre of `C`'s hyperelliptic map; `C` has a rational
Weierstrass point, so `C : y^2 = q1(u)` with `q1` a quintic and that fibre at `u = 0`, and `Y` is
`y^2 = q1(t^2/c)` for a twist `c` supported on the bad primes. The star's equation is the one fitted
to its differentials (`star-model.log`, 58 good primes against the newform), a candidate: nothing
here proves it, and the genus-four model inherits that status. Run from the repository root:

| script | what it shows |
|---|---|
| `enumerate.m` | over every fibre of the star above a rational point of height ≤ 10^6 and every twist supported on {−1, 2, 5, 7, 29}, which double covers have the point counts of `Y`: exactly one, the fibre `x = 1/2` with `c = 1` |
| `thuemahler.m` | the branch fibre, exactly: the two branch points stay distinct modulo every prime of good reduction, so the weighted y-coordinate of the fibre is a unit outside {2, 5, 7, 29} and the quintic factor of the star's sextic is an S-unit, a Thue–Mahler equation whose 17 solutions (`thuemahler.sols`, Gherga–Siksek's solver, 246 s) give 6 possible fibres, of which only `x = 1/2` passes the count test |
| `verify.m` | the resulting curve, its two genus-two quotients (the star and `C'`), its counts against the trace formula at the 31 good primes below 150, and the conductor of `C'` |
