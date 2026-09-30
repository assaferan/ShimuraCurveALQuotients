// Chunked enumeration for odd q: all y^2 = f, deg f in {2g+1,2g+2}, squarefree, constant term f(0) = c0 = CHUNK.
// Writes the set of Weil polys (char poly form) to OUT; union over chunks compared with the table afterwards.
SetColumns(0);
g := StringToInteger(g); q := StringToInteger(q); c0 := StringToInteger(CHUNK);
F := GF(q); P<x> := PolynomialRing(F);
S := {};
key := func<W | "[" cat Join([IntegerToString(c) : c in Reverse(Coefficients(W))], ",") cat "]">;
for v in CartesianPower(F, 2*g+2) do
  f := P!([F!c0] cat [v[i] : i in [1..2*g+2]]);
  if Degree(f) lt 2*g+1 or Discriminant(f) eq 0 then continue; end if;
  Include(~S, key(Reverse(LPolynomial(HyperellipticCurve(f)))));
end for;
Fo := Open(OUT, "w"); for s in S do fprintf Fo, "%o\n", s; end for; delete Fo;
quit;
