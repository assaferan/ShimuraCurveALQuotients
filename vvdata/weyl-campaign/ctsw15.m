// GAP 5, step 1: is a_E(m) a class number times LOCAL FACTORS?
//
// E is now an explicit mass-weighted combination of ternary genus thetas
// (Thm 9.11), so by Siegel-Weil each term is an Eisenstein series whose
// coefficients are (archimedean) x prod_p beta_p(m).  The prediction is
//
//     a_E(m) = H(4m) x prod_{p | DN} lambda_p(m)
//
// with lambda_p depending only on ord_p(m) and a quadratic symbol -- the exact
// analogue at the m-th coefficient of what c_p^ram and (c_p^Eich - c_p^ram)/2
// are at the constant term.  This tests it at 15_2.
//
// Magma has no HurwitzClassNumber, so H is built here and GATED against the
// standard table before use.
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/ctlat.m";
DEPTH := 200;

// Hurwitz class number H(n), n >= 0
Hur := function(n)
    if n eq 0 then return Rationals()!(-1)/12; end if;
    if n mod 4 in {1, 2} then return Rationals()!0; end if;
    s := Rationals()!0;
    f := 1;
    while f^2 le n do
        if n mod f^2 eq 0 then
            d := -(n div f^2);
            if d mod 4 in {0, 1} then
                w := d eq -3 select 6 else (d eq -4 select 4 else 2);
                s +:= Rationals()!(2*ClassNumber(QuadraticForms(d)))/w;
            end if;
        end if;
        f +:= 1;
    end while;
    return s;
end function;

// GATE: Gauss's three-square identity, with r_3 computed rather than recalled.
//   r_3(n) = 12 H(4n)   for n = 1,2 mod 4
//   r_3(n) = 24 H(n)    for n = 3 mod 8
// (A hardcoded H-table was tried first and the gate fired on H(32) -- the table
// was wrong, not the routine.  Gate against an identity, not against memory.)
T3 := ThetaSeries(StandardLattice(3), 400);
r3 := func< n | Coefficient(T3, n) >;
bad := [ ];
for n in [1..90] do
    if n mod 4 in {1,2} and r3(n) ne 12*Hur(4*n) then Append(~bad, <n,"12H(4n)">); end if;
    if n mod 8 eq 3 and r3(n) ne 24*Hur(n) then Append(~bad, <n,"24H(n)">); end if;
end for;
printf "Hurwitz gate (Gauss three-square identity): %o failures%o\n", #bad,
    #bad eq 0 select "  -- routine trusted" else Sprintf(" %o", bad);
error if #bad ne 0, "Hurwitz class number routine is wrong -- stop";

thetabar := function(Dp, Rl)
    A := ChangeRing(grossGram(Dp, Rl), Rationals());
    L := LatticeWithGram(A);
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

// 15_2 : D = 15 = 3*5, N = 2, s = 3, D' = 5
v1 := thetabar(5, 2); vs := thetabar(5, 6); v3 := thetabar(30, 1);
G := [ Rationals() | -(1/2)*v1[i] + vs[i] - (1/2)*v3[i] : i in [1..DEPTH] ];

printf "\n  m   a_E(m)    H(4m)     ratio      (-m|3) (-m|5)  v2 v3 v5\n";
ratios := AssociativeArray();
for m in [1..DEPTH-1] do
    if G[m+1] eq 0 then continue; end if;
    h := Hur(4*m);
    r := h eq 0 select Rationals()!0 else G[m+1]/h;
    k3 := KroneckerSymbol(-m, 3); k5 := KroneckerSymbol(-m, 5);
    key := <k3, k5, Valuation(m,2), Valuation(m,3), Valuation(m,5)>;
    if not IsDefined(ratios, key) then ratios[key] := {@ @}; end if;
    Include(~ratios[key], r);
    if m le 60 then
        printf "%3o  %7o  %8o  %9o     %2o     %2o    %o  %o  %o\n",
            m, G[m+1], h, r, k3, k5,
            Valuation(m,2), Valuation(m,3), Valuation(m,5);
    end if;
end for;

printf "\n=== is a_E(m)/H(4m) constant on each local class? ===\n";
nconst := 0; ntot := 0;
for key in Sort([ k : k in Keys(ratios) ]) do
    ntot +:= 1;
    if #ratios[key] eq 1 then nconst +:= 1; end if;
    printf "  (-m|3)=%2o (-m|5)=%2o  v2=%o v3=%o v5=%o :  %o distinct  %o\n",
        key[1], key[2], key[3], key[4], key[5], #ratios[key],
        #ratios[key] le 4 select [ x : x in ratios[key] ] else "(many)";
end for;
printf "\n%o of %o local classes carry a CONSTANT ratio\n", nconst, ntot;
printf "DONE\n";
quit;
