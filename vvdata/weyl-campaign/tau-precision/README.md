# tau-precision: why M0MultiplierExact loses accuracy as M grows

Probes from the 2026-09-02/03 investigation of the `10_61` gate-3 failure. Full write-up:
`PLAN.md` on `tier1-models` (REPAIR track).

All of these need **no Borcherds forms and no pipeline** -- `VVCosetReps` / `VVSTWord` /
`slashdata` and `ds = Divisors(M)` depend only on `M = 2DN` (`4DN` odd). ~7 min at `M = 1220`.

    tauprobe.m    per-coset Im(z0) at one base; ranks the failing cosets
    tauscale.m    the M^2 law across four bases
    tauwindow.m   eta dynamic range per coset + the safe-window band table
    amtest.m      unrelated: the scalar a_E vs measured A_m comparison (1 of 13)

## What they found

`min Im(z0) = Im(tau0)/M^2 * (1 + O(1/M))`, verified at M = 60/580/748/1220 with the ratio
converging to 1 (1.0011 at M = 1220). The repo's own recorded precision -- exact, 33 digits,
18 digits, catastrophic -- tracks that collapse in lockstep, so "achievable precision is
base-dependent" has a cause: `M`, through the coset denominators.

Safe window (worst eta spread per Im(z0) band, all 2232 cosets at M = 1220):

    [0, 1e-5) 1137 -> 92.17    [1e-5,1e-4) 669 -> 37.45    [1e-4,1e-3) 307 -> 9.90
    [1e-3,1e-2)  91 ->  3.55   [1e-2,1e-1)  22 ->  9.25    [1e-1,1e0)    5 -> 100.26
    [1e0, 1e2)    1 -> 181.56

Both ends fatal, middle safe. `wi = 2` fails at the LARGE end: `d` runs to `M`, so
`Im(d*z0) = 882` and eta underflows to 1e-100 -- predicted -100.31, measured -100.272.

## CAVEAT -- read before quoting any small-Im number

`tauwindow.m`'s `logeta` is a **one-S-step approximation**: it uses `-C*Im(z)` when `Im(z) >= 0.05`
and otherwise one application of `eta(-1/z) = sqrt(-iz) eta(z)`. That is valid only when `Im(z)` is
large, or when a single `S` makes it large. When BOTH `z` and `-1/z` sit near the real axis --
i.e. `z` near a rational with large denominator -- neither branch applies and the value is wrong.

`wi = 1221` (`Im(z0) = 8.8e-7`) is exactly that regime, so its reported 1.39-digit spread is
**provisional**. Re-measure with a real `DedekindEta` before drawing any conclusion about it.
