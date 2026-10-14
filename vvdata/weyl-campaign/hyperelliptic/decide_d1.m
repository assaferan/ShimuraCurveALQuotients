// D = 1: the holomorphic differentials of X_0(N)/W are the weight-2 cusp forms on Gamma_0(N) fixed
// by the Atkin-Lehner involutions in W, so Noether's count runs on q-expansions: the number of
// linear relations among the pairwise products of a basis is (g-1)(g-2)/2 for a hyperelliptic
// curve and (g-2)(g-3)/2 otherwise.  Precision past the Sturm bound of weight 4 on Gamma_0(N).
SetColumns(0);
todo := [<84,[12,21,28],3>, <96,[32],3>, <96,[96],3>, <114,[2,19,38],3>, <120,[3,5,15],3>, <130,[2,13,26],3>, <130,[5,26,130],3>, <156,[12,39,52],3>, <165,[5,11,55],4>, <180,[5,36,180],5>, <180,[20,36,45],5>, <195,[5,39,195],3>, <210,[3,10,14,30,35,42,105],3>, <210,[2,5,7,10,14,35,70],3>, <210,[5,6,7,30,35,42,210],3>, <240,[3,5,15,16,48,80,240],3>, <294,[2,3,6,49,98,147,294],3>, <330,[6,10,11,15,66,110,165],4>, <510,[2,3,5,6,10,15,17,30,34,51,85,102,170,255,510],3>, <570,[2,3,5,6,10,15,19,30,38,57,95,114,190,285,570],4>, <630,[2,5,7,9,10,14,18,35,45,63,70,90,126,315,630],5>, <714,[2,3,6,7,14,17,21,34,42,51,102,119,238,357,714],5>];
for t in todo do
    N := t[1]; W := t[2]; g := t[3];
    S := CuspForms(Gamma0(N), 2);
    V := VectorSpace(Rationals(), Dimension(S));
    fixed := V;
    for q in W do
        Wq := AtkinLehnerOperator(S, q);
        fixed := fixed meet Kernel(Wq - 1);
    end for;
    dimfix := Dimension(fixed);
    // Sturm bound for weight 4
    idx := N * &*[1 + 1/p : p in PrimeDivisors(N)];
    prec := Integers()!(4 * idx / 12) + 10;
    B := Basis(S);
    qexp := [PowerSeries(b, prec) : b in B];
    R := Parent(qexp[1]);
    fs := [&+[v[i] * qexp[i] : i in [1..#B]] : v in Basis(fixed)];
    prods := [fs[i] * fs[j] : j in [i..#fs], i in [1..#fs]];
    M := Matrix(Rationals(), [[Coefficient(f, n) : n in [0..prec-1]] : f in prods]);
    nrel := #prods - Rank(M);
    hyp := nrel eq (g-1)*(g-2) div 2; non := nrel eq (g-2)*(g-3) div 2;
    printf "X_0(%o)/<w_%o>: genus %o (fixed space dimension %o), %o quadrics: %o\n", N, W, g, dimfix, nrel,
        hyp select "HYPERELLIPTIC" else (non select "not hyperelliptic" else "NEITHER COUNT");
end for;
quit;
