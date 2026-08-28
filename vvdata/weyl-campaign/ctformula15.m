// IS E* LITERALLY THE GROSS COMBINATION?
//
// The a0_g route needs a HOLOMORPHIC weight-3/2 form whose constant term at
// every cusp is (M/g) rho_g; at 15_2 it was FOUND by linear algebra and
// rationalised to the eta quotient
//     E_eis = 4/5 mono1 - 4 mono2 - 4/5 mono3 + 4 mono4      (eisrat, ce4f18e)
// The N>1 theorem says the mass-weighted Gross combination
//     -1/2 thetabar(5,2) + thetabar(5,6) - 1/2 thetabar(30,1)
// has exactly those constant terms -- and genus theta series ARE holomorphic
// weight-3/2 forms.  So the theorem CONSTRUCTS E*, rather than asserting it
// exists, at every base in its scope.
//
// The two agree in constant terms by the theorem.  This script asks the sharper
// question: do they agree as Q-EXPANSIONS?  If yes, E* is the Gross combination
// on the nose and W(m) is a combination of ternary representation numbers.  If
// they differ by a cusp form, W(m) is unaffected (Borcherds obstruction: f's
// principal part pairs to zero with every weight-3/2 cusp form, because f
// exists) -- but then the difference is itself worth naming.
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/ctlat.m";
DEPTH := 60;
R<q> := PowerSeriesRing(Rationals(), DEPTH);

eta_unit := function(d0, e)
    ser := R!1; n := 1;
    while d0*n lt DEPTH do ser *:= (1 - q^(d0*n))^e; n +:= 1; end while;
    return ser;
end function;
ds := [1,2,3,4,5,6,10,12,15,20,30,60];
mono := function(rr)
    s1 := &+[ ds[i]*rr[i] : i in [1..#ds] ];
    ser := R!1;
    for i in [1..#ds] do
        if rr[i] ne 0 then ser *:= eta_unit(ds[i], rr[i]); end if;
    end for;
    return q^(s1 div 24) * ser;
end function;
E := (4/5)*mono([-6,15,0,-6,0,0,0,0,0,0,0,0])
     - 4*mono([-3,7,0,-3,1,0,0,0,0,1,0,0])
     - (4/5)*mono([-1,0,0,4,-1,0,3,0,0,-2,0,0])
     + 4*mono([-1,2,0,0,0,0,0,0,-3,0,7,-2]);

// mass-weighted average theta of the Gross genus, normalised to constant term 1
thetabar := function(Dp, Rl)
    G := ChangeRing(grossGram(Dp, Rl), Rationals());
    L := LatticeWithGram(G);
    reps := Representatives(Genus(L));
    raw := [ Rationals() | 0 : m in [0..2*DEPTH] ];
    den := Rationals()!0;
    for Lr in reps do
        w := Rationals()!1/#AutomorphismGroup(Lr);
        T := ThetaSeries(Lr, 2*DEPTH);
        for m in [0..2*DEPTH] do raw[m+1] +:= w*Coefficient(T, m); end for;
        den +:= w;
    end for;
    raw := [ x/den : x in raw ];
    // detect the norm convention: if every odd index vanishes, the quadratic
    // form is half the Gram norm and the series lives on even indices
    halved := forall{ m : m in [1..2*DEPTH by 2] | raw[m+1] eq 0 };
    printf "   thetabar(%o,%o): %o classes in genus, det %o, %o\n",
        Dp, Rl, #reps, Determinant(G),
        halved select "odd indices vanish -> reading q^(m/2)" else "all indices live";
    if halved then
        return [ raw[2*m+1] : m in [0..DEPTH-1] ];
    end if;
    return [ raw[m+1] : m in [0..DEPTH-1] ];
end function;

printf "=== 15_2: is E* the Gross combination? ===\n";
v1 := thetabar(5, 2);
vs := thetabar(5, 6);
v3 := thetabar(30, 1);
G := [ Rationals() | -(1/2)*v1[i] + vs[i] - (1/2)*v3[i] : i in [1..DEPTH] ];
Ec := [ Coefficient(E, m) : m in [0..DEPTH-1] ];

printf "\n m :   E*(eta quotient)      Gross combination      difference\n";
for m in [0..Minimum(29, DEPTH-1)] do
    printf " %2o : %14o %22o %18o\n", m, Ec[m+1], G[m+1], Ec[m+1]-G[m+1];
end for;
d := [ Ec[i] - G[i] : i in [1..DEPTH] ];
printf "\ndifference: %o\n",
    forall{ x : x in d | x eq 0 } select "IDENTICALLY ZERO -- E* IS the Gross combination"
    else Sprintf("nonzero, first at m = %o, support %o",
        Minimum([ i-1 : i in [1..DEPTH] | d[i] ne 0 ]),
        [ i-1 : i in [1..DEPTH] | d[i] ne 0 ]);
printf "DONE\n";
quit;
