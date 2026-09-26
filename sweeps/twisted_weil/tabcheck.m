// Independent enumeration of ALL hyperelliptic curves of genus g over F_q (q = p), compare the set of
// Weil polynomials with data/hypg<g>q<q>.txt.  Char p odd: y^2 = f, deg f in {2g+1, 2g+2}, f squarefree.
// Char 2: y^2 + h y = f with deg h <= g+1, deg f <= 2g+2, smooth genus g (HyperellipticCurve + Genus).
SetColumns(0);
g := StringToInteger(g); q := StringToInteger(q);
F := GF(q); P<x> := PolynomialRing(F);
S := {};
key := func<W | "[" cat Join([IntegerToString(c) : c in Reverse(Coefficients(W))], ",") cat "]">;
if q ne 2 then
  for v in CartesianPower(F, 2*g+3) do
    f := P![v[i] : i in [1..2*g+3]];
    if Degree(f) lt 2*g+1 or Discriminant(f) eq 0 then continue; end if;
    C := HyperellipticCurve(f);
    Include(~S, key(Reverse(LPolynomial(C))));   // char poly t^{2g} L(1/t)
  end for;
else
  for vh in CartesianPower(F, g+2) do
    h := P![vh[i] : i in [1..g+2]];
    if h eq 0 then continue; end if;
    for vf in CartesianPower(F, 2*g+3) do
      f := P![vf[i] : i in [1..2*g+3]];
      ok := IsHyperellipticCurve([f, h]);
      if not ok then continue; end if;
      C := HyperellipticCurve(f, h);
      if Genus(C) ne g then continue; end if;
      Include(~S, key(Reverse(LPolynomial(C))));
    end for;
  end for;
end if;
T := Set([l : l in Split(Read(Sprintf("data/hypg%oq%o.txt", g, q)), "\n") | #l gt 0]);
printf "g=%o q=%o enumerated %o table %o equal %o  enum-minus-table %o table-minus-enum %o\n", g, q, #S, #T, S eq T, #(S diff T), #(T diff S);
quit;
