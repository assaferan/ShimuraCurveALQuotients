// Noether's count over a list of quotients, as a pipeline stage -- OFF unless asked for.
//
// ⚠⚠ WHY THIS STAGE DOES NOT SET `IsSubhyp`, AND MUST NOT.  Every other filter in the pipeline
// marks a curve with a verdict that is PROVED: the twisted Weil polynomials are exact, the dual
// graph is combinatorial, the trace-formula counts are exact.  What this stage computes is the
// numerical kernel dimension of a matrix of quadratic monomials in floating-point values of
// differentials, and calling a singular value zero is not a proof that it is zero.  Writing such
// a reading into the field that the rest of the pipeline treats as settled would put a number of
// one kind where every consumer expects another -- the failure this project has paid for more
// than once.  So the stage touches no curve record at all: it returns its readings, and the
// caller decides what to do with evidence.
//
// WHAT IT IS FOR.  The readings predict what the exact routes should find, which is worth a great
// deal when an exact search is expensive: at X_0^46(31) the reading said hyperelliptic, and a
// widened divisor search later proved it, unconditionally.  A reading that disagrees with an
// exact verdict is a defect in one of them and should be chased, not averaged.
//
// COST, and why it is conditional.  A verdict is Points evaluations of Theorem 3.1, each a sum
// over the lattice points of a quaternary form, and the cost is linear in D*N: seconds at D*N in
// the tens, about 700 s per evaluation at D*N = 2274 and 90 digits.  A sweep is CPU-days, so the
// stage is inert unless Enabled is set or NELSON is in the environment.
//
// READING THE RESULT.  At genus g the canonical image lies on (g-1)(g-2)/2 quadrics when the
// curve is hyperelliptic and (g-2)(g-3)/2 when it is not: one against none at genus 3, three
// against one at genus 4.  A verdict is accepted only when it is stable, and stability means
// insensitivity to the PRECISION, not to the cutoff -- two cutoffs at one precision can share a
// precision floor and agree to many digits while both are wrong (X_0^4130(1), which reads
// non-hyperelliptic at 50 digits at three cutoffs and hyperelliptic at 80).  Because a floor puts
// a nonzero number where a vanishing one belongs, the reading it fakes is the NEGATIVE one, so a
// "NOT hyperelliptic" answer is the one that needs the second precision; NelsonHyperellipticVerdict
// enforces that through ConfirmPrecision.

intrinsic NelsonFilterReadings(curves::SeqEnum[ShimuraQuot] :
                               Cutoff := 40, Precision := 80, Points := 18,
                               ConfirmPrecision := 0, MinGenus := 3, MaxLevel := 0) -> SeqEnum
{Noether's count on each quotient of genus at least MinGenus, returned as tuples
<CurveID, D, N, W, genus, verdict, smallest singular value, gap, relations, lattice terms>.  Quotients of
smaller genus are skipped, as are those with D*N above MaxLevel when that is nonzero, which is how
a caller bounds the cost.  Nothing is written to the curve records: see the header for why.}
    out := [];
    for X in curves do
        if not assigned X`g or X`g lt MinGenus then continue; end if;
        if MaxLevel gt 0 and X`D * X`N gt MaxLevel then continue; end if;
        vprintf ShimuraQuotients, 1 : "Nelson count: X_0^%o(%o)/%o, genus %o\n",
                X`D, X`N, Sort(Setseq(X`W)), X`g;
        ok := true; verdict := "ERROR"; rec := 0;
        try
            verdict, rec := NelsonHyperellipticVerdict(X`D, X`N, X`W :
                                Cutoff := Cutoff, Precision := Precision, Points := Points,
                                ConfirmPrecision := ConfirmPrecision);
        catch e
            ok := false;
            vprintf ShimuraQuotients, 1 : "  refused: %o\n", e`Object;
        end try;
        if ok then
            Append(~out, <X`CurveID, X`D, X`N, Sort(Setseq(X`W)), X`g, verdict,
                          Minimum(rec`singular), rec`gap, rec`relations, rec`terms>);
        else
            Append(~out, <X`CurveID, X`D, X`N, Sort(Setseq(X`W)), X`g, "REFUSED", 0, 0, -1, 0>);
        end if;
    end for;
    return out;
end intrinsic;

intrinsic FilterByNelsonCount(~curves::SeqEnum : Enabled := false, Cutoff := 40, Precision := 80,
                              Points := 18, ConfirmPrecision := 0, MaxLevel := 0, Out := "")
{The pipeline hook for Noether's count.  INERT unless Enabled is true or NELSON is set in the
environment, because a sweep costs CPU-days.  It marks no curve: the readings are numerical and
the pipeline's verdict fields are for proved ones, so they are printed and, with Out, written to
that file for a human to weigh against the exact routes.}
    if not (Enabled or GetEnv("NELSON") ne "") then
        vprintf ShimuraQuotients, 2 : "Nelson count: not enabled, skipping\n";
        return;
    end if;
    rs := NelsonFilterReadings(curves : Cutoff := Cutoff, Precision := Precision, Points := Points,
                                        ConfirmPrecision := ConfirmPrecision, MaxLevel := MaxLevel);
    for r in rs do
        vprintf ShimuraQuotients, 1 : "  X_0^%o(%o)/%o  genus %o  %o  (%o relations, smallest %o, gap %o)\n",
                r[2], r[3], r[4], r[5], r[6], r[9], r[7], r[8];
    end for;
    if Out ne "" then
        lines := ["// Noether's count, NUMERICAL readings -- not proofs, and not verdicts of the",
                  "// pipeline.  Columns: CurveID, D, N, W, genus, reading, smallest singular value,",
                  "// gap, relations, lattice terms.  See the header of NelsonFilter.m."];
        for r in rs do
            Append(~lines, Sprintf("<%o, %o, %o, %o, %o, \"%o\", %o, %o, %o, %o>,",
                   r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9], r[10]));
        end for;
        Write(Out, Join(lines, "\n") : Overwrite);
    end if;
end intrinsic;
