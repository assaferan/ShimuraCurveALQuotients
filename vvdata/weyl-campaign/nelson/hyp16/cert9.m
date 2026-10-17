AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
D := 39; N := 11; DN := 429;
NFs := Newforms(CuspForms(Gamma0(DN), 2));
printf "all newform orbits of level exactly %o, weight 2:\n", DN;
printf "%o\n", "deg | a_3,a_11,a_13 | quaternionic w_3,w_11,w_13,w_39,w_429";
tot := 0;
for orb in NFs do
    f := orb[1];
    if Level(Parent(f)) ne DN then continue; end if;
    dg := Degree(BaseRing(Parent(f)));
    q := qExpansion(f, 14);
    a3 := Coefficient(q,3); a11 := Coefficient(q,11); a13 := Coefficient(q,13);
    // quaternionic: a_p at p | D (=3,13), -a_p at p | N (=11)
    w3 := a3; w13 := a13; w11 := -a11; w39 := a3*a13; w429 := a3*a13*(-a11);
    inW := (w11 eq 1) and (w39 eq 1) and (w429 eq 1);
    printf "%-6o %-24o %-28o W-fixed=%-6o (contributes %o)\n", dg, [a3,a11,a13], [w3,w11,w13,w39,w429], inW, inW select dg else 0;
    if inW then tot +:= dg; end if;
end for;
printf "\ntotal dimension of the W-fixed space: %o   (genus of X/W is 3)\n", tot;
printf "\nthe two rational classes of conductor 429:\n";
for F in EllipticCurves(CremonaDatabase(), 429) do
    r := CremonaReference(F);
    if r ne "429a1" and r ne "429b1" then continue; end if;
    printf "   %-7o a_3=%-3o a_11=%-3o a_13=%-3o\n", r, TraceOfFrobenius(F,3), TraceOfFrobenius(F,11), TraceOfFrobenius(F,13);
end for;
print "CERT9 DONE";
