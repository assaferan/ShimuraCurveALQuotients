# The 122 never-started bases, triaged

Predictor-only sweep — `ppint.m` (gate 1, principal-part integrality) then `cmsupply.m` (gate 2,
CM supply). No model generation. Run on **lovelace** (256 cores, idle) at 64 concurrent with a
3600 s cap per tool. Full per-base table and logs alongside; tally in `SUMMARY.txt`.

     81  CAPPED-1h      no verdict even at 3600 s
     21  OBSTRUCTED     Borcherds obstruction
     13  VX-ASSERT      the `assert vx ge 0` method boundary
      3  ASSERT
      2  NONINTEGRAL
      1  CM-STARVED     206_3, "Could not find enough rational CM points"
      1  INTEGRAL
    122  TOTAL

## The two results that matter

**Runnable candidates — 2 of 122.**

    10_61   INTEGRAL     CM OK margin 0     run as-is
    14_43   NONINTEGRAL  CM OK margin 0     needs IntegralSolution (the fixable gate)

Both sit at **margin 0** on CM supply — exactly enough points, no slack. That is the position
`34_11` was in, and it cleared gates 1-3 before stopping at class-constancy.

**The obstructed class grows 28 -> 49.** Twenty-one more bases fail with "Failed to find all
Borcherds forms":

    106_5 10_59 118_5 122_5 14_41 14_47 202_3 214_3 218_3 226_3 22_29
    254_3 26_23 34_19 38_17 46_13 58_11 6_103 82_7 86_7 94_7

Per [[borcherds-obstruction-is-real]] neither more divisor triples nor deeper poles can help any
of them, and per [[even-correction-blocked-on-theory]] the escape hatch needs a theorem that does
not yet exist. **That theorem is now worth 49 bases, not 28** — which materially changes its
priority relative to everything else in the backlog.

## The methodological finding: these predictors are NOT cheap here

`ppint` was minutes-per-base on the NONINTEGRAL class. On this cohort **81 of 122 produced no
verdict in a full hour**. Those bases were never started precisely because they are large, and
`ppint` must construct Borcherds forms before it can say anything about integrality — so on this
cohort it is an hours-scale tool, not a triage tool.

A first pass at a 600 s cap was abandoned for this reason: it converted "slow" into "unknown"
faster without buying resolution. **The cap was bounding the measurement itself.** Raising it to
3600 s converted 60 unknowns into 41 real verdicts, so the longer cap was worth it, but the
remaining 81 need a fundamentally cheaper predictor rather than more wall-clock.

## Traps hit while running this (all mine)

* **macOS has no `timeout`**, so the local pass capped with `kill` — and Magma BUFFERS stdout to
  a file, so a killed run leaves a **0-byte log**, indistinguishable from "not started". Eleven
  of twelve local bases produced nothing and were reported as "running". GNU `timeout` on
  lovelace truncates cleanly instead; the remote pass records an explicit `SWEEPCAP` marker.
* **`pkill -x magma` kills nothing** — the binary is `magma.exe` behind a wrapper. Worse, the
  verification used the same broken pattern and so reported success. A 48-way pass launched on
  top of a still-running 32-way pass, both writing to one directory with different caps; all 54
  logs from that window were quarantined rather than trusted. **Verify a kill with a different
  pattern than the one used to kill.**
* Committed cache files carry escaped `\[` line-starts, so a byte-comparison against freshly
  generated Normaliz output fails even when the vector sets are identical. Compare as sets.
