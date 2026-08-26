# SCOPING: what X_0^1155(1) costs, and in what order to attack it

Companion to `prereg-1155.md` (the predictions). This file is the
engineering brief: the three walls, their measured sizes, and the
design for the one that needs a rewrite. Nothing here has been run
against the base; everything is either measured on existing bases or
computed from the lattice.

Base facts, all measured:

    D = 1155 = 3*5*7*11,  N = 1,  M(eta) = 2DN? NO -- odd DN, so M = 4DN = 4620
    |A| = 2,668,050 = 2 * 3^2 * 5^2 * 7^2 * 11^2
    Smith invariants of the Gram: [1, 1155, 2310]
    2-part: Z/2 of ODD type  =>  the r2=1 branch of eis32k applies unchanged
    #divisors(4620) = 48
    #cosets of Gamma_0(4620) = 13,824

Comparison base is `330_1`: `M = 660`, `|A| = 217,800`, 24 divisors,
1728 cosets, 171-member pool.

## WALL A -- the eta pool at M = 4620 (do this first; it gates everything)

No banked base has `M = 4620`, so unlike the odd-level quartet (which
reused pools by level) this one must be enumerated.

**Solved already: the character key.** `enum32m` needs
`(s1 mod 24, s2 mod 24, parity bits per prime)`, taken from a MONO
dump. Measured across ALL six banked pools it is universal:

    s1 = s2 = 0 mod 24,  parity bit 1 at p = 2 only

(verified on epool_15_2/21_2/33_2/35_2/210_1/330_1, i.e. both parities
of `DN` -- the odd-level quartet reused these pools and fitted, so the
key transfers to odd `DN`). For `M = 4620` the primes are
`[2,3,5,7,11]`, so the key is `(0, 0, (1,0,0,0,0))`. Build
`mono1155_1_synth.log` on the model of `mono330_1_synth.log`: the BASE
line plus one carrier MONO vector satisfying that key.

**The cost.** `enum32m` is meet-in-the-middle over supports of size
`<= K = 7 = 3 + 4` with entries in `[-R, R] = [-8, 8]` (`nv = 17`):

| | A-cache rows `C(nd,3)*nv^3` | B-loop rows `C(nd,4)*nv^4` |
|---|---|---|
| `M = 660` (nd = 24) | 9.94e6 | 8.88e8 |
| `M = 4620` (nd = 48) | 8.50e7 (**8.5x**) | 1.63e10 (**18.3x**) |

so roughly **18x the 660 pool time** and **8.5x its peak memory** --
the A-cache alone is order 2 GB of int64 at `nd = 48`, which is the
likelier failure mode.

**Partial calibration, measured.** `enum32m.py mono330_1_synth.log <out>
8 7` at `M = 660` ran **>25 minutes without finishing** and was killed
(no output file). So `R=8, K=7` at `M = 4620` is **>7 hours** on the
same machine, and plausibly far worse once the A-cache starts swapping.
Caveat: the parameters the banked `epool_330_1.txt` was actually
generated with are not recorded, so this may be a harsher setting than
was ever used in production -- re-derive from a completed 660 run
before trusting the multiplier.

**Conclusion: do not attempt `R=8, K=7` at 4620.** Start from the
reduced settings below and only widen if the pool fails to span.

**If it is too big**, the cheap levers in order: drop `R` from 8 to 6
(`nv` 17 -> 13, B-loop x0.34); drop `K` from 7 to 6 (`NA=3, NB=3`,
B-loop `C(48,3)*17^3` = 8.5e7, a 190x cut); prune the A-cache to
subsets that can still reach `sum r = 3`. A smaller pool is fine as
long as it spans -- the residual is the arbiter, and `330_1` fitted
with 171 members against 1728 cosets.

## WALL B -- the O(|A|) setup in eis32k (the rewrite, and it is worth doing)

At `|A| = 2.67e6` (12.25x `330_1`) these all enumerate the whole group:

* `VVWeilFFT` builds `elts`
* `lift` -- `|A|` rational 3-vectors via `@@Ld\`to_disc`
* `Qbar` -- `|A|` rationals, one quadratic-form evaluation each
* `idx` -- an AssociativeArray with `|A|` keys (the memory hazard)
* `Syl[p]` -- a full scan per Jordan prime (5 of them here)
* `YTab` -- `O(|A|)` per key, ~23 keys at 330_1
* the `a = 0` word -- one `VVRhoInvE0FFT` over `|A|`

**The rewrite is easy and the payoff is large, because the components
are tiny.** The Jordan decomposition here is

    J_2 = Z/2 (2)   J_3 = 9   J_5 = 25   J_7 = 49   J_11 = 121

-- the largest is 121 elements. Everything the driver needs is an
orthogonal sum over these, so nothing ever has to touch 2.67e6 objects:

* `Qbar(x) = sum_p Qbar_p(x_p) mod 1` and
  `B(x,y) = sum_p B_p(x_p, y_p)` -- the Jordan decomposition is
  orthogonal (already relied on by `gsumJ`, which is per-component).
* `|A[c]| = prod_p |J_p[c]| = prod_p |J_p[p^v_p(c)]|` -- replaces the
  `dcn` scan.
* `c*y = tgt` splits componentwise: in `J_p` write `c = p^v * u` with
  `u` a `p`-adic unit; solvable iff `tgt_p` lies in `p^v J_p`, and then
  `y_p = u^{-1} (tgt_p / p^v)`. This replaces the `y0` search.
* the tracked coset `est` and `x_c` are specified componentwise (`x_c`
  lives in `J_2` alone).
* `gsumJ` needs no change at all.

So represent a coset as a tuple of component indices, and drop `elts`,
`lift`, `idx` entirely. This also speeds every existing base, and the
`CTL:=1` gate validates it at 15_1/15_2/22_3 in seconds before it is
trusted at 1155.

**The `a = 0` word is the one piece left.** It is currently measured by
a single FFT because that was "one word, affordable"; at 2.67e6 that
one call may dominate the whole run. Either try it and time it, or
derive the entry in closed form -- for the bare-`S` coset
`(rho(S)^{-1} e_0)_eta` is a Gauss-sum normalisation, independent of
`eta`, so the only real content is the sign the code currently measures.

## WALL C -- the E-table (likely the actual wall-clock sink)

`eis32k` evaluates `a0at(w, rE, 3)` for every pool member at every
coset. `330_1` was already E-table-dominant (see
`eis32k-closed-form-driver`). At 4620 the cosets go 1728 -> **13,824**
(8x) and the pool will be at least as large, so budget **>= 8-16x the
330_1 E-table**. It is embarrassingly parallel over pool members;
splitting the pool across processes and concatenating the EMAT dump is
the obvious move, and `kernelrat`/`eisrat` consume the dump either way.

## SUGGESTED ORDER

1. **Wall B rewrite first**, gated by `CTL:=1` on 15_1 / 15_2 / 22_3.
   It is self-contained, it does not need the pool, and it speeds every
   future base. Do NOT start it at 1155 -- validate small.
2. **Calibrate Wall A** by regenerating the 660 pool and timing it;
   then decide `R`/`K` for 4620 and build `mono1155_1_synth.log`.
3. Build the 4620 pool.
4. Run `eis32k DD:=1155 NN:=1 PR:=120 EF:=<pool>` -- parallelise the
   E-table if step 3's pool is large.
5. `kernelrat` -> `eisrat` -> `genallgross` -> `allgross1155_1`.
6. **`genrestrict` on all five candidate supports** (Q4 is the
   falsifier and needs every one of them, not just the winner):
   `385,1;385,3` / `105,1;105,11` / `11,1;11,105` /
   `231,1;231,5` / `165,1;165,7`.
7. Scorecard against `prereg-1155.md`.

## GOTCHAS THAT WILL BITE

* `PR:=120` at `M >= 660`; the `a0at` consistency check is relative,
  not absolute, for exactly this reason.
* `mpmath` and `numpy` are needed by kernelrat/eisrat/genallgross and
  enum32m; scratchpad venvs purge -- recreate.
* Run magma from the campaign ROOT (AttachSpec), paths as
  `vvdata/weyl-campaign/...`.
* zsh: `${=var}` to word-split; quote `echo "==="`; `timeout` is not
  installed on this machine.
* Flatten the 80-column magma wrap before parsing any log.
* `M = 2DN` is the eta level; the weight-3/2 form level is `4DN`. At
  odd `DN` the eta level IS `4DN`. Conflating them already cost one
  preregistered prediction.
