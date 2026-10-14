// Controls for the q-expansion count: Ogg's list of hyperelliptic X_0(N) (genus >= 3: 40, 48 of
// genus 3; 47 of genus 4; 59 of genus 5) against non-hyperelliptic ones of the same genera (34, 45,
// 64 of genus 3; 53 of genus 4; 55 of genus 5), all with W trivial.
SetColumns(0);
function count(N, W)
    S := CuspForms(Gamma0(N), 2);
    V := VectorSpace(Rationals(), Dimension(S)); fixed := V;
    for q in W do fixed := fixed meet Kernel(AtkinLehnerOperator(S, q) - 1); end for;
    g := Dimension(fixed);
    idx := N * &*[1 + 1/p : p in PrimeDivisors(N)];
    prec := Integers()!(4 * idx / 12) + 10;
    B := Basis(S); qexp := [PowerSeries(b, prec) : b in B];
    fs := [&+[v[i] * qexp[i] : i in [1..#B]] : v in Basis(fixed)];
    prods := [fs[i] * fs[j] : j in [i..#fs], i in [1..#fs]];
    nrel := #prods - Rank(Matrix(Rationals(), [[Coefficient(f, n) : n in [0..prec-1]] : f in prods]));
    return g, nrel, nrel eq (g-1)*(g-2) div 2;
end function;
for N in [40, 48, 47, 59, 34, 45, 64, 53, 55] do
    g, nrel, hyp := count(N, []);
    printf "X_0(%o): genus %o, %o quadrics, hyperelliptic %o (Ogg: %o)\n", N, g, nrel, hyp, N in {40, 48, 47, 59};
end for;
// a quotient with a known answer: X_0(114)/<w_2, w_38> has an explicit genus-3 model y^2 = f(x) (PR #51)
g, nrel, hyp := count(114, [2, 38]); printf "X_0(114)/<w_2,w_38>: genus %o, %o quadrics, hyperelliptic %o (expected true)\n", g, nrel, hyp;
quit;
