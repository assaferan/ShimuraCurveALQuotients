// Which Zagier pullbacks Z_d live in M_{3/2}(120) at all, and is E_eis in
// Gross + (those Z_d) + cusp with new relations appearing?
DEPTH := 200; MMAX := DEPTH - 1;
lev := 120; DN := 30;
M := HalfIntegralWeightForms(lev, 3/2);
BM := Basis(M, DEPTH);
AM := Matrix(Rationals(), MMAX+1, #BM, [ [Coefficient(BM[j], m) : j in [1..#BM]] : m in [0..MMAX] ]);
rkM := Rank(AM);
function H24(k)
    if k eq 0 then return Rationals()!(-1/12); end if;
    if k mod 4 in {1,2} then return Rationals()!0; end if;
    tot := Rationals()!0;
    for f in [1..Isqrt(k)] do
        if k mod (f*f) ne 0 then continue; end if;
        d := -(k div (f*f));
        if ((d mod 4) + 4) mod 4 notin {0,1} then continue; end if;
        w := d eq -3 select 1/3 else (d eq -4 select 1/2 else 1);
        tot +:= w*ClassNumber(d);
    end for;
    return tot;
end function;
for d in Divisors(DN) do
    col := Matrix(Rationals(), MMAX+1, 1, [ (m mod d eq 0) select H24(m div d) else 0 : m in [0..MMAX] ]);
    printf "Z_%o in M_{3/2}(120): %o\n", d, Rank(HorizontalJoin(AM, col)) eq rkM;
end for;
quit;
