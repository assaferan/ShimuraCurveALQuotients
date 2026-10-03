# Weil-stage re-timing, 2026-10-02

`FilterByWeilPolynomial` replayed prime by prime on **lovelace with the class-number tables**, on
`main` at `95b19e6` (so with #56, #57 and #58), over **36 curves** spanning genus 3 to 7, `#W` from
1 to 64 and level from 30 to 30030 -- including the five curves the 2026-09-30 run had timed.

    weil_design.m      picks the curves (one per genus / #W / level-bucket cell) and prints each
                       one's admissible primes and trace-formula term count
    weil_retime.sh     the launcher used on lovelace (8 at a time, 12 h limit each)
    weil_timing.m      the driver, unchanged from the 2026-09-30 run
    logs/              per-prime CPU and wall times, one file per curve
    retime.psv         id | D,N,g,#W,Qmax | total cpu s | primes done, scraped from the logs
    fit.py, fit.log    the comparison and the power-law fit

## What it says

⚠ **Read the warning at the end first**: an earlier version of this file was written from 34 of the
36 curves, and the two missing ones were the two slowest of the set.

**The stage reaches about 1.7 h on a single curve.** The two slowest are `1416` = `X_0(595)/W8`
(`g = 5`) at 101 min and `2325` = `X_0^6(97)/W2` (`g = 6`) at 74 min; the next are 33 and 32 min.
`1071` takes 6 min, against 87 min with #58 alone and 14.8 h before it, so #56 and #57 are another
factor of 14 on that curve.

**The cost is `#W * sum_p p^(g/2)`**, the shape of the dominant trace at `n = p^g`:
Eichler-Selberg sums over the `w` in `W`, and each inner sum runs over `t^2 < 4n/Q_w`. (The
actual term count carries `Q_w^(-1/2)` weights, and that version fits worse, so this is the shape
of the sum rather than its length.) Of the 271 curve pairs whose times differ by more than a factor 10 it
orders 246 correctly, against 198 for the sum of `p^g` that `CurveCostProxy` used before; at a
factor 2, 440 of 520 against 343.

**That is as much as the data supports.** Among `#W * sum p^(a g)` for `a` from 0.5 to 0.8, and
`#W * max p^(a g)`, the inversion counts differ by less than the noise of 36 curves (63 to 97 of
520 at the factor-2 threshold, non-monotone in `a`). The term count is kept because it is the
principled form, not because it scores best. Dropping `#W` is clearly worse (192 of 520), and so is
the exponent `g` in place of `g/2` (177 of 520).

**⚠ The heavy shape is high genus at a MODERATE level — not a large `W`.** The two slowest curves
have `#W = 8` and `#W = 2`, both at level about 590. `#W` earns its place in the ordering, not in
the extremes. The level's fitted exponent is `+0.15` (`#W^a (sum p^(g/2))^b (D N)^c`), so weakly
positive rather than absent.

**⚠ A residual spread of 237x remains and real pairs are mis-ordered** -- `2325` (74 min) sits
below `13029` (86 s). This is an ordering heuristic for the parallel dispatch, never a cost. Quote
the measured times.

## ⚠ What the first version of this file got wrong

It was written when `1416` and `2325` had not finished, and they were the two slowest curves of the
set. On 34 curves the estimate appeared to order 216 of 220 pairs with a spread of 34x, the level's
exponent came out `-0.07`, and the heavy shape looked like high genus with a LARGE `W`. All four
statements changed once the two finished. The lesson is the one in `CLAUDE.md` about truncated runs:
a batch that is 34/36 done looks finished, and here the stragglers were the whole point.
