SetColumns(0);
DN := 429; ps := [2,5,7,17,19,23];
printf "rational newform orbits of level exactly 429:\n";
for orb in Newforms(CuspForms(Gamma0(DN), 2)) do
    f := orb[1];
    if Level(Parent(f)) ne DN or Degree(BaseRing(Parent(f))) ne 1 then continue; end if;
    q := qExpansion(f, 32);
    printf "   a_2,a_3,a_5,a_7,a_11,a_13 = %o\n",
        [ Integers() ! Coefficient(q,n) : n in [2,3,5,7,11,13] ];
end for;
printf "\nelliptic curves, one per isogeny class of conductor 429:\n";
seen := {};
for F in EllipticCurves(CremonaDatabase(), DN) do
    r := CremonaReference(F); lbl := r[1..4];
    if lbl in seen then continue; end if; Include(~seen, lbl);
    printf "   %-6o a_2,a_3,a_5,a_7,a_11,a_13 = %o\n", lbl,
        [ TraceOfFrobenius(F,n) : n in [2,3,5,7,11,13] ];
end for;
E := MinimalModel(EllipticCurve([0,0,0,-24,63]));
printf "\nour genus-1 quotient's Jacobian: %o\n", CremonaReference(E);
printf "   a_2,a_3,a_5,a_7,a_11,a_13 = %o\n", [ TraceOfFrobenius(E,n) : n in [2,3,5,7,11,13] ];
print "CERT11 DONE";
