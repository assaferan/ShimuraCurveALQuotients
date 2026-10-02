# Weil-stage re-timing, 2026-10-02

`FilterByWeilPolynomial` replayed prime by prime on **lovelace with the class-number tables**, on
`main` at `95b19e6` (so with #56, #57 and #58), over **34 curves** spanning genus 3 to 7, `#W` from
1 to 64 and level from 30 to 30030 -- including the five curves the 2026-09-30 run had timed.

    weil_design.m      picks the curves (one per genus / #W / level-bucket cell) and prints each
                       one's admissible primes and trace-formula term count
    weil_retime.sh     the launcher used on lovelace (8 at a time, 12 h limit each)
    weil_timing.m      the driver, unchanged from the 2026-09-30 run
    logs/              per-prime CPU and wall times, one file per curve
    retime.psv         id | D,N,g,#W,Qmax | total cpu s | primes done, scraped from the logs
    fit.py, fit.log    the comparison and the power-law fit

## What it says

**The stage is minutes per curve now, not hours.** The slowest of the 34 took 33 min. The three
curves of the earlier run came in at 9.5 min (`7296`), 6 min (`1071`) and 86 s (`13029`); `1071`
had taken 14.8 h on `main` before #58 and 87 min with #58 alone, so #56 and #57 account for
another factor of 14.

**The cost is `#W * sum_p p^(g/2)`.** That is the number of terms the dominant trace at `n = p^g`
needs: Eichler-Selberg sums over the `w` in `W`, and each inner sum runs over `t^2 < 4n/Q_w`, about
`2 sqrt(n/Q_w)` values. Of the 220 curve pairs whose times differ by more than a factor 10, it
orders 216 correctly; the sum of `p^g`, which `CurveCostProxy` used before, orders 166.

**The level does not enter.** Fitting `#W^a (sum p^(g/2))^b (D N)^c` gives `c = -0.07`: no
measurable effect over levels from 30 to 30030. So `Qmax` belongs in the prime bounds, not in the
cost.

**⚠ The `Q_w^(-1/2)` weights make the fit WORSE.** The exact term count
`sum_p sum_w 2 sqrt(p^g / Q_w)` mis-orders 90 of the 452 pairs differing by more than 2x, against
44 for the unweighted `#W * sum p^(g/2)`. So the cost grows with `#W` rather than shrinking with
the `Q_w`: the per-term work is not constant, and a larger `W` costs more than its shorter `t`-ranges
save. The best simple fit is `#W^1.3 (sum p^(g/2))^1.6`, which cuts the residual spread from 34x to
24x -- not used in the proxy, which only has to order the chunks.

**⚠ A residual spread of 34x remains**, so this is an ordering heuristic and not a cost model. Quote
the measured times, never the estimate.

Two curves of the 36 launched (`1416`, `2325`) had not finished when the fit was taken; both are
mid-range and neither changes the comparison.
