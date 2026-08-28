// GAP A, step 1: is the Gross combination IN the holomorphic eta pool at 15_2?
//
// If it is, then G is an eta quotient combination like E*, the whole eis32b
// machinery applies to it unchanged, and A_f = sum_w a0((f G)|w) can be
// computed with the existing code.  If it is not, the eta span is too small and
// the comparison has to be done vector-valued.
//
// Also re-derives E* in the same basis as a control: the banked eisrat answer
// is 4/5 mono1 - 4 mono2 - 4/5 mono3 + 4 mono4, so the solver must find a
// solution and the difference of the two solutions must be the cusp form
// measured in ctformula15.log.
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/ctlat.m";
DEPTH := 60;
R<q> := PowerSeriesRing(Rationals(), DEPTH);
ds := Divisors(60);
printf "ds = %o\n", ds;

eta_unit := function(d0, e)
    ser := R!1; n := 1;
    while d0*n lt DEPTH do ser *:= (1 - q^(d0*n))^e; n +:= 1; end while;
    return ser;
end function;
mono := function(rr)
    s1 := &+[ ds[i]*rr[i] : i in [1..#ds] ];
    error if s1 mod 24 ne 0, "monomial not integral in q", rr;
    ser := R!1;
    for i in [1..#ds] do
        if rr[i] ne 0 then ser *:= eta_unit(ds[i], rr[i]); end if;
    end for;
    return q^(s1 div 24) * ser;
end function;

Epool := [ ];
for r in eval Read("vvdata/weyl-campaign/epool_15_2.txt") do
    Append(~Epool, [ Integers() | x : x in r ]);
end for;
printf "#pool = %o\n", #Epool;
P := [ mono(r) : r in Epool ];

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
    return [ raw[2*m+1] : m in [0..DEPTH-1] ];
end function;

v1 := thetabar(5, 2); vs := thetabar(5, 6); v3 := thetabar(30, 1);
Gv := [ Rationals() | -(1/2)*v1[i] + vs[i] - (1/2)*v3[i] : i in [1..DEPTH] ];
Ev := [ Coefficient((4/5)*mono([-6,15,0,-6,0,0,0,0,0,0,0,0])
                    - 4*mono([-3,7,0,-3,1,0,0,0,0,1,0,0])
                    - (4/5)*mono([-1,0,0,4,-1,0,3,0,0,-2,0,0])
                    + 4*mono([-1,2,0,0,0,0,0,0,-3,0,7,-2]), m) : m in [0..DEPTH-1] ];

A := Matrix(Rationals(), DEPTH, #P, [ [ Coefficient(P[j], m) : j in [1..#P] ] : m in [0..DEPTH-1] ]);
for nm->tgt in [* <"E* (eisrat)", Ev>, <"Gross combination", Gv> *] do
    lab := tgt[1]; v := Vector(Rationals(), tgt[2]);
    ok, sol := IsConsistent(Transpose(A), v);
    printf "\n%-20o in the eta pool span: %o\n", lab, ok;
    if ok then
        printf "   coefficients: %o\n", [ sol[j] : j in [1..#P] | sol[j] ne 0 ];
        printf "   on pool rows: %o\n", [ j : j in [1..#P] | sol[j] ne 0 ];
        chk := [ &+[ sol[j]*Coefficient(P[j], m) : j in [1..#P] ] : m in [0..DEPTH-1] ];
        printf "   reconstruction exact: %o\n", chk eq tgt[2];
    end if;
end for;
printf "pool span rank = %o of %o\n", Rank(A), #P;
printf "DONE\n";
quit;
