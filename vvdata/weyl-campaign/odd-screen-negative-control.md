# The even screen is anti-correlated with the truth on odd `D` — 6 controls

Measured 2026-09-14, read off `~/shimura/defic/ODD_*.log`, recorded 2026-09-15.

The even-`D` screen decides on `deficit = Ncols - Rank`, swept over the pole order `P`. Applying
that same statistic to six **odd**-`D` bases that all build and all have committed models:

    base   ladder over P                        verdict the even rule would give
    15_1   2  4  7 11 15 20   (P = 10 ... 266)  OBSTRUCTED  -- wrong
    21_2   2  3  4  5  8 10   (P = 17 ... 266)  OBSTRUCTED  -- wrong
    39_1   4  5  7 10 14 19   (P = 37 ... 266)  OBSTRUCTED  -- wrong
    51_1   6  9  9 15 20      (P = 51 ... 266)  OBSTRUCTED  -- wrong
    55_1   7  8 10 16 22      (P = 52 ... 266)  OBSTRUCTED  -- wrong
    57_1   5  7 10 15 19      (P = 53 ... 266)  OBSTRUCTED  -- wrong

**6 of 6 false.** The point is not that the even statistic is noisy on odd `D`. It is that the
deficit *rises monotonically with `P`* on every one of these, so more evidence makes the wrong
verdict look stronger — the opposite of the even case, where the diagnostic is invariance across
rungs. Any rule of the form "large invariant deficit ⇒ obstructed" is therefore not merely
unreliable on odd `D`, it is actively misleading.

This is why the odd ladder is over `m` and the statistic is `wdef` (the deficit restricted to the
achievable targets) — see `deficit_odd.m` and `HANDOFF.md`, 2026-09-14 (later). `55_1` reads
`deficit 3 / wdef 0` and builds; `33_1` and `69_1` read plain deficit 1 and 3 and both build.

⚠ The converse limit still holds and is not touched by this control: an odd `wdef >= 1` is **not**
a verdict either (`21_2`, a Guo-Yang base with a committed model, exhausts its whole `m` ladder at
`wdef >= 2`). **Odd screens give CLEAR verdicts only.**
