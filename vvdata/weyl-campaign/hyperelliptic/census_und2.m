// The undecided curves at bases without a model file: level M of the Borcherds search, its number
// of divisors (the Normaliz wall is at 24), the genus of the star curve (the method needs 0), and
// whether the base has any genus-0 quotient at all.
AttachSpec("ShimuraQuotients.spec"); SetColumns(0);
curves := eval Read("data/curves_after_UpdateCurves8.dat");
und := [X : X in curves | not assigned X`IsSubhyp];
bases := {<X`D, X`N> : X in und | not FileExists(Sprintf("data/models/models_%o_%o.m", X`D, X`N))};
rows := [];
for b in Sort(SetToSequence(bases)) do
    D := b[1]; N := b[2]; DN := D*N;
    D0 := DN; while IsEven(D0) do D0 div:= 2; end while;
    M := 4*D0;
    full := {Integers()| d : d in Divisors(DN) | GCD(d, DN div d) eq 1};
    star := [Y : Y in curves | Y`D eq D and Y`N eq N and Y`W eq full];
    gstar := #star eq 1 select star[1]`g else -1;
    n0 := #[Y : Y in curves | Y`D eq D and Y`N eq N and Y`g eq 0];
    nu := #[X : X in und | X`D eq D and X`N eq N];
    Append(~rows, <D, N, M, #Divisors(M), gstar, n0, nu>);
end for;
printf "bases without a file: %o, undecided curves there: %o\n", #rows, &+[r[7] : r in rows];
printf "star curve of genus 0: %o bases (%o curves); star genus > 0: %o bases (%o curves)\n",
    #[r : r in rows | r[5] eq 0], &+[Integers()| r[7] : r in rows | r[5] eq 0], #[r : r in rows | r[5] gt 0], &+[Integers()| r[7] : r in rows | r[5] gt 0];
printf "star genus 0 and #div(M) < 24 (below the wall): %o bases, %o curves; #div(M) >= 24: %o bases, %o curves\n",
    #[r : r in rows | r[5] eq 0 and r[4] lt 24], &+[Integers()| r[7] : r in rows | r[5] eq 0 and r[4] lt 24],
    #[r : r in rows | r[5] eq 0 and r[4] ge 24], &+[Integers()| r[7] : r in rows | r[5] eq 0 and r[4] ge 24];
printf "below the wall, star genus 0 (D, N, M, #div M, g*, genus-0 quotients, undecided): %o\n", [r : r in rows | r[5] eq 0 and r[4] lt 24];
printf "star genus > 0 (D, N, M, #div M, g*, genus-0 quotients, undecided): %o\n", [r : r in rows | r[5] gt 0];
printf "on the wall (D, N, M, #div M, g*, genus-0 quotients, undecided): %o\n", [r : r in rows | r[5] eq 0 and r[4] ge 24];
quit;
