// Read hypsearch candidates (C c a0..a8), keep squarefree f, print the set of Weil polys (charpoly form).
SetColumns(0); q := StringToInteger(q); F := GF(q); P<x> := PolynomialRing(F);
key := func<W | "[" cat Join([IntegerToString(c) : c in Reverse(Coefficients(W))], ",") cat "]">;
S := {* *}; ns := 0;
for l in Split(Read(IN), "\n") do f := Split(l, " "); if #f lt 3 or f[1] ne "C" then continue; end if;
  c := F!StringToInteger(f[2]); g := c * P![F!StringToInteger(f[i]) : i in [3..#f]];
  if Discriminant(g) eq 0 then ns +:= 1; continue; end if;
  Include(~S, key(Reverse(LPolynomial(HyperellipticCurve(g)))));
end for;
printf "singular %o smooth-by-poly %o\n", ns, S; quit;
