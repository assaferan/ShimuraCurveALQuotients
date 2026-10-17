AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
DN := 429;  // D = 39 = 3*13, N = 11; quaternionic w_p = a_p at p|D, -a_p at p|N
plusdim := 0; minusdim := 0; plusforms := [];
for orb in Newforms(CuspForms(Gamma0(DN), 2)) do
    f := orb[1];
    if Level(Parent(f)) ne DN then continue; end if;
    dg := Degree(BaseRing(Parent(f)));
    q := qExpansion(f, 14);
    a3 := Coefficient(q,3); a11 := Coefficient(q,11); a13 := Coefficient(q,13);
    w3 := a3; w11 := -a11; w39 := a3*a13; w429 := a3*a13*(-a11);
    if not (w11 eq 1 and w39 eq 1 and w429 eq 1) then continue; end if;
    printf "W-fixed orbit of degree %o : a_3 = %o, a_11 = %o, a_13 = %o, so w_3 = %o\n", dg, a3, a11, a13, w3;
    if w3 eq 1 then plusdim +:= dg; Append(~plusforms, f); else minusdim +:= dg; end if;
end for;
printf "\nW-fixed space: %o + %o = %o  (genus of X/W is 3)\n", plusdim, minusdim, plusdim+minusdim;
printf "w_3 = +1 part has dimension %o, so X/<W,w_3> has genus %o -- the model says 1\n", plusdim, plusdim;
if #plusforms eq 1 then
    f := plusforms[1]; q := qExpansion(f, 30);
    E := MinimalModel(EllipticCurve([0,0,0,-24,63]));
    ps := [2,5,7,17,19,23,29];
    af := [ Integers() ! Coefficient(q,p) : p in ps ];
    aE := [ TraceOfFrobenius(E,p) : p in ps ];
    printf "\nthe w_3 = +1 form   a_p for p in %o: %o\n", ps, af;
    printf "429b1, our Jacobian a_p for p in %o: %o\n", ps, aE;
    printf "IDENTICAL: %o\n", af eq aE;
end if;
print "CERT10 DONE";
