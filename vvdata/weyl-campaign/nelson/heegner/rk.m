SetColumns(0);
DB := CremonaDatabase();
for cond in [1326, 2730] do
    printf "\nconductor %o:\n", cond;
    for i in [1..NumberOfIsogenyClasses(DB, cond)] do
        E := EllipticCurve(DB, cond, i, 1);
        ap := func<p | TraceOfFrobenius(E, p)>;
        tw := cond eq 1326 select [8, -3, -51, 17, -136] else [8, -15, -39, -4, -105, -273];
        printf "  %o: a_p at %o = %o  rank %o tors %o | twist ranks %o\n", CremonaReference(E), PrimeDivisors(cond), [ap(p) : p in PrimeDivisors(cond)], Rank(E), #TorsionSubgroup(E), [<d, Rank(QuadraticTwist(E, d))> : d in tw];
    end for;
end for;
