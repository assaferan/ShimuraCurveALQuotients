// Which arithmetic sequences lie in M_{3/2}(Gamma_0(60)) (Magma's convention)?
// Membership battery: rank test of each candidate against the space's basis.
DEPTH := 120;
MMAX := DEPTH - 1;
R<q> := PowerSeriesRing(Rationals(), DEPTH);

M := HalfIntegralWeightForms(60, 3/2);
BM := Basis(M, DEPTH);
nb := #BM;
AM := Matrix(Rationals(), MMAX+1, nb, [ [Coefficient(BM[j], m) : j in [1..nb]] : m in [0..MMAX] ]);
rk := Rank(AM);
printf "dim M = %o (rank %o to depth %o)\n", nb, rk, DEPTH;

// Hurwitz table
HMAX := 225*MMAX;
HTAB := [Rationals()| 0 : i in [1..HMAX+1]];
for N in [1..HMAX] do
    if N mod 4 in {1,2} then continue; end if;
    tot := 0;
    for f in [f : f in [1..Isqrt(N)] | N mod f^2 eq 0] do
        d := -(N div f^2);
        if d mod 4 in {0,1} then
            h := ClassNumber(QuadraticForms(d));
            ww := d eq -3 select 6 else (d eq -4 select 4 else 2);
            tot +:= 2*h/ww;
        end if;
    end for;
    HTAB[N+1] := tot;
end for;
function HH(n) return n lt 0 select 0 else HTAB[n+1]; end function;

// theta series
function theta(a)
    ser := R!0;
    n := 0;
    while a*n^2 lt DEPTH do ser +:= (n eq 0 select 1 else 2)*q^(a*n^2); n +:= 1; end while;
    return ser;
end function;

cands := [* *]; names := [];
// H-shifts (with H(0) = -1/12 convention for the constant term)
for pair in [<1,1>,<4,1>,<1,3>,<4,3>,<1,5>,<4,5>,<1,15>,<4,15>,<1,2>,<1,4>,<1,6>,<1,10>,<1,12>,<1,20>,<1,30>,<1,60>,<4,2>,<4,6>,<4,10>,<4,30>] do
    mult := pair[1]; d := pair[2];
    seq := [Rationals()| (mult*m mod d eq 0) select HH(mult*m div d) else 0 : m in [0..MMAX] ];
    seq[1] := -1/12 * (d eq 1 select 1 else 0);  // constant term of the Zagier series under V_d/U
    Append(~cands, seq); Append(~names, Sprintf("H(%om/%o)", mult, d));
end for;
// U_{p^2}-images: H(9m/d), H(25m/d), H(225m/d)
for pair in [<9,1>,<9,3>,<25,1>,<25,5>,<225,1>,<225,15>] do
    mult := pair[1]; d := pair[2];
    seq := [Rationals()| (mult*m mod d eq 0) select HH(mult*m div d) else 0 : m in [0..MMAX] ];
    seq[1] := -1/12 * 0;
    Append(~cands, seq); Append(~names, Sprintf("H(%om/%o)", mult, d));
end for;
// ternary diagonal thetas with a,b,c in divisors(15)
divs15 := [1,3,5,15];
for a in divs15 do for b in divs15 do for c in divs15 do
    if a gt b or b gt c then continue; end if;
    ser := theta(a)*theta(b)*theta(c);
    Append(~cands, [Rationals()| Coefficient(ser, m) : m in [0..MMAX] ]); Append(~names, Sprintf("th%o.th%o.th%o", a, b, c));
end for; end for; end for;

for i in [1..#cands] do
    bcol := Matrix(Rationals(), MMAX+1, 1, cands[i]);
    inM := Rank(HorizontalJoin(AM, bcol)) eq rk;
    printf "%-16o : %o\n", names[i], inM select "IN" else "out";
end for;
quit;
