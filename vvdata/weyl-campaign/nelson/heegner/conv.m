// Which sign convention for the Atkin-Lehner eigenvalues of E = X/W' holds: E's eigenvalues must be +1 on all of
// W' for some index-2 overgroup W' of W with genus 1. Convention A: eps_v = a_v for v | D, -a_v for v | N.
// Convention B: eps_v = -a_v everywhere. Convention C: eps_v = a_v everywhere.
AttachSpec("ShimuraQuotients.spec"); SetColumns(0);
load "/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/heeg/rowsW.m";
star := func<a, b | a*b div GCD(a, b)^2>;
for c in cases do
    D, N, W, lab := Explode(c); W := {Integers() | w : w in W};
    E := EllipticCurve(lab);
    if Conductor(E) ne D*N then continue; end if;
    als := [d : d in Divisors(D*N) | GCD(d, (D*N) div d) eq 1 and d notin W];
    overs := {}; for u in als do Wu := W join {star(u, w) : w in W}; if GenusShimuraCurveQuotient(D, N, Wu) eq 1 then Include(~overs, Wu); end if; end for;
    res := "";
    for conv in ["A", "B", "C"] do
        eps := func<v | conv eq "A" select ((D mod v eq 0) select TraceOfFrobenius(E, v) else -TraceOfFrobenius(E, v)) else (conv eq "B" select -TraceOfFrobenius(E, v) else TraceOfFrobenius(E, v))>;
        epsm := func<m | &*[Integers() | eps(v) : v in PrimeDivisors(m)]>;
        good := [Wu : Wu in overs | &and[epsm(m) eq 1 : m in Wu]];
        res cat:= Sprintf(" %o:%o", conv, #good gt 0 select Setseq(good[1]) else "none");
    end for;
    printf "%o_%o W=%o E=%o a_p=%o |%o\n", D, N, Setseq(W), lab, [TraceOfFrobenius(E, p) : p in PrimeDivisors(D*N)], res;
end for;
