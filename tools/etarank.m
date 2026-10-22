// Rank of the q-expansions at oo of a set of eta quotients, independently of the repository's code:
// the spanning checks of the shift sets in polymake/README.md.
//   magma -b M:=780 file:=<exponent vectors> npole:=472 tools/etarank.m < /dev/null
// For a derived file (first line a comment) and for a polytope file alike.
// Each row of the file is an exponent vector r indexed by the divisors d of M (ascending);
// the eta quotient is prod_d eta(d tau)^{r_d} = q^{(sum d r_d)/24} prod_d prod_n (1 - q^{d n})^{r_d},
// and the polytope's congruence sum d r_d = 0 mod 24 makes the exponent of q an integer, so pole
// orders are in units of q.  Rank over GF(p) is a lower bound for the rank over Q, and equals it
// unless p divides a minor; we use two primes.  Precision: a weight-1/2 form on Gamma_0(M),
// holomorphic away from oo, with a pole of order at most npole at oo has at most
// index/24 + npole zeros at oo, so coefficients up to q^(index/24 + npole + 1) determine the form.
M := StringToInteger(M); npole := StringToInteger(npole);
pts := eval Read(file);
ds := Divisors(M);
assert #ds eq #pts[1];
index := M * &*[1 + 1/p : p in PrimeDivisors(M)];
T := Ceiling(index/24) + npole + 2;
printf "M=%o  file=%o  %o vectors  precision q^%o\n", M, file, #pts, T;
for p in [10007, 1000003] do
    F := GF(p);
    R<q> := PowerSeriesRing(F, T + npole + 1);
    E := AssociativeArray();
    for d in ds do
        e := R!1;
        n := 1;
        while d*n le T + npole do e *:= (1 - q^(d*n)); n +:= 1; end while;
        E[d] := e;
    end for;
    pw := AssociativeArray();
    rows := [];
    leading := [];
    for r in pts do
        v24 := &+[ds[i]*r[i] : i in [1..#ds]];
        assert v24 mod 24 eq 0;
        v := v24 div 24;                              // order at oo in units of q
        assert v ge -npole;
        f := R!1;
        for i in [1..#ds] do
            if r[i] ne 0 then
                if not IsDefined(pw, <i, r[i]>) then pw[<i, r[i]>] := E[ds[i]]^r[i]; end if;
                f *:= pw[<i, r[i]>];
            end if;
        end for;
        f := q^(v + npole) * f;                       // shift so every row starts at q^0
        Append(~rows, [Coefficient(f, j) : j in [0..T + npole]]);
        Append(~leading, v);
    end for;
    A := Matrix(F, rows);
    Ech := EchelonForm(A);
    rk := Rank(A);
    piv := [Min([j : j in [1..Ncols(Ech)] | Ech[i,j] ne 0]) - 1 - npole : i in [1..rk]];
    gaps := [v : v in [-npole..Max(piv)] | v notin piv];
    printf "p=%o: rank %o; leading orders from %o to %o; %o orders below the top not attained (gaps), largest gap at %o\n",
           p, rk, Min(piv), Max(piv), #gaps, #gaps gt 0 select Max(gaps) else "none";
end for;
