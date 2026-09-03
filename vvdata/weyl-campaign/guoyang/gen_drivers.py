#!/usr/bin/env python3
"""Generate one self-contained Magma verification driver per (D,N) base, checking the
PRIMARY hauptmodule ('s', the full-quotient star curve) column of its Guo-Yang CM-value
table via a cross-ratio (Mobius-invariant) comparison -- no assumption that our hauptmodule
shares Guo-Yang's zero/pole normalisation, only that 3 reference points are shared.
"""
import json

with open("tables.json") as f:
    tables = json.load(f)

OUT = "drivers"
import os
os.makedirs(OUT, exist_ok=True)

TEMPLATE = '''SetQuitOnError(true);
AttachSpec("ShimuraQuotients.spec");
D := {D}; N := {N};

// (disc, GY published value) pairs -- rational or Infinity() -- Guo-Yang "CM-values of
// X_0^{D}({N})", primary hauptmodule column (full star quotient).
gy := [
{rows}
];

function mobius(z0, z1, z2, z)
    num := (z eq Infinity() or z0 eq Infinity()) select 1 else z - z0;
    den := (z eq Infinity() or z1 eq Infinity()) select 1 else z - z1;
    c1  := (z2 eq Infinity() or z1 eq Infinity()) select 1 else z2 - z1;
    c2  := (z2 eq Infinity() or z0 eq Infinity()) select 1 else z2 - z0;
    if den*c2 eq 0 then return Infinity(); end if;
    return (num*c1)/(den*c2);
end function;

t0 := Realtime();
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){{c : c in curves | IsStarCurve(c)}};

keep := {{ t[1] : t in gy }};
tab := ValuesAtCMPoints(star, curves : Keep := keep);
discs := tab`Discs;
srow := tab`Values[tab`sIndex];
idx := AssociativeArray();
for i->d in discs do idx[d] := i; end for;

printf "VERIFY D=%o N=%o : %o gy rows, table has %o discs\\n", D, N, #gy, #discs;
for t in gy do
    d := t[1];
    if not IsDefined(idx, d) then
        printf "  MISSING disc %o (Keep did not admit it)\\n", d;
    end if;
end for;

// reference triple: first 3 rows with distinct GY values and present in our table
have := [t : t in gy | IsDefined(idx, t[1])];
error if #have lt 4, "fewer than 4 usable rows";
ref := [];
for t in have do
    if #ref eq 3 then break; end if;
    ok := true;
    for r in ref do if r[2] eq t[2] then ok := false; end if; end for;
    if ok then Append(~ref, t); end if;
end for;
error if #ref lt 3, "could not find 3 rows with distinct GY values";
z0 := srow[idx[ref[1][1]]]; z1 := srow[idx[ref[2][1]]]; z2 := srow[idx[ref[3][1]]];
w0 := ref[1][2]; w1 := ref[2][2]; w2 := ref[3][2];
printf "  reference discs %o %o %o\\n", ref[1][1], ref[2][1], ref[3][1];

nchecked := 0; nfail := 0;
for t in have do
    d, expected := Explode(t);
    if d in {{ref[1][1], ref[2][1], ref[3][1]}} then continue; end if;
    got := mobius(z0, z1, z2, srow[idx[d]]);
    want := mobius(w0, w1, w2, expected);
    if got ne want then
        printf "  FAIL disc %o: got %o, want %o (published %o)\\n", d, got, want, expected;
        nfail +:= 1;
    end if;
    nchecked +:= 1;
end for;
printf "RESULT D=%o N=%o checked=%o failed=%o time=%os\\n", D, N, nchecked, nfail, Realtime()-t0;
quit;
'''

count = 0
for t in tables:
    D, N = t["D"], t["N"]
    # A discriminant can legitimately label MORE THAN ONE CM point on the curve (distinct
    # embeddings with the same discriminant), and Guo-Yang's tables then list it twice with
    # different values (seen at 10_19, disc -52: rows infinity and -4). Our lookup is by
    # discriminant alone and can't disambiguate which point a repeated disc means, so such
    # discriminants are excluded rather than risk a spurious mismatch.
    from collections import Counter
    disc_counts = Counter(r["disc"] for r in t["rows"])
    dup_discs = {d for d, n in disc_counts.items() if n > 1}
    rows = []
    for r in t["rows"]:
        if r["disc"] in dup_discs:
            continue
        v = r["vals"][0]
        if v["kind"] == "rat":
            num, den = v["val"]
            if den == 1:
                rows.append(f"<{r['disc']}, {num}>")
            else:
                rows.append(f"<{r['disc']}, {num}/{den}>")
        elif v["kind"] == "oo":
            rows.append(f"<{r['disc']}, Infinity()>")
    if len(rows) < 4:
        continue
    content = TEMPLATE.format(D=D, N=N, rows=",\n".join(rows))
    fname = f"{OUT}/verify_{D}_{N}.m"
    with open(fname, "w") as f:
        f.write(content)
    count += 1

print(f"wrote {count} driver scripts to {OUT}/")
