// X_0^1155(1): the unrestricted Gross-span test AND all five preregistered
// candidate supports, in ONE process.
//
// TWO FIXES over the generated drivers.
//
// (1) DEPTH.  allgross_template.m hard-codes DEPTH := 200, which fixes the
// number of ROWS of the rank-test matrix.  Its COLUMNS are the GROSS averages
// plus the ENTIRE cuspidal basis, and at level 4620 dim S = 561, so the
// unrestricted matrix is 200 x 601 and each restricted one is 200 x 563 --
// more columns than rows.  The cusp columns alone then reach rank 200 and
// EVERY candidate "passes" for free: the test cannot distinguish anything.
// That would have looked like the prereg's most interesting outcome ("several
// supports carry E_eis, so the support is a convention") while being pure
// artifact.  DEPTH = 800 > 601 restores content, and the guard below makes the
// requirement explicit instead of assumed.
//   This never bit the earlier bases -- dim S is 85 at 210_1 and 133 at 330_1,
//   both under 200 -- and the banked verdicts are DEPTH-STABLE: re-running
//   330_1 at DEPTH = 800 reproduces rank 148 vs 148, aliasing dim 25, E in
//   span: true, exactly as banked.  So this is a fix for THIS base, not a
//   correction to the campaign.
//
// (2) ONE construction.  Basis(M, DEPTH) cost is in CONSTRUCTING the basis,
// not in the depth -- measured at level 1320: Basis(200) 84.6 s, Basis(400)
// 0.37 s, Basis(800) 11.8 s, memory flat.  Running the six tests as six
// drivers would repeat that construction six times.  Here the space, its
// basis, and all 40 GROSS columns are built ONCE, and the six rank tests are
// column selections off them.
DEPTH := 800;   // see the guard below: must EXCEED ngross + dim S
MMAX := DEPTH - 1;
R<q> := PowerSeriesRing(Rationals(), DEPTH);
DN := 1155;
dsv := [1,2,3,4,5,6,7,10,11,12,14,15,20,21,22,28,30,33,35,42,44,55,60,66,70,77,84,105,110,132,140,154,165,210,220,231,308,330,385,420,462,660,770,924,1155,1540,2310,4620];
rs := [ [-2,5,0,-2,0,0,-4,0,0,0,10,0,0,0,0,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [-2,5,0,-2,0,0,0,0,-4,0,0,0,0,0,10,0,0,0,0,0,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [-1,2,0,0,0,0,-3,0,0,0,7,0,0,0,0,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [-1,2,0,0,0,0,-1,0,0,0,0,0,0,2,0,1,0,0,0,1,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [-1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,1,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0],
        [-1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-3,0,0,0,0,7,0,0,-2,0,0,0,0],
        [0,-1,0,2,0,0,0,0,-2,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,-1,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,1,0,0,0,3,0,-1,0],
        [0,0,-2,0,0,5,0,0,-3,-2,0,0,0,0,6,0,0,1,0,0,0,0,0,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,-2,0,0,0,5,0,0,0,2,-2,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,2,0,0,-1,-1,0,-1,0,0,0,0,0,0,0,0],
        [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,2,0,2,-1,0,0,-1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0] ];
cs := [ -1/10, 1/9, -3/5, 8/5, -48/5, 48/5, -8/9, -12, 2/9, 3/5, 3/2, -4/3 ];
lev := 4*DN;
printf "BASE 1155_1: DN = %o, level %o\n", DN, lev;

function eta_unit(d0, e)
    ser := R!1; n := 1;
    while d0*n lt DEPTH do ser *:= (1 - q^(d0*n))^e; n +:= 1; end while;
    return ser;
end function;
function mono(rr, dsx)
    s1 := &+[dsx[i]*rr[i] : i in [1..#dsx]];
    error if s1 mod 24 ne 0, "bad monomial";
    ser := R!1;
    for i in [1..#dsx] do
        if rr[i] ne 0 then ser *:= eta_unit(dsx[i], rr[i]); end if;
    end for;
    return q^(s1 div 24) * ser;
end function;
E := &+[ cs[i]*mono(rs[i], dsv) : i in [1..#rs] ];

M := HalfIntegralWeightForms(lev, 3/2);
S := CuspidalSubspace(M);
printf "dim M = %o, dim S = %o (Eis %o)\n", Dimension(M), Dimension(S), Dimension(M)-Dimension(S);
BM := Basis(M, DEPTH);
BS := Basis(S, DEPTH);
AM := Matrix(Rationals(), MMAX+1, #BM, [ [Coefficient(BM[j], m) : j in [1..#BM]] : m in [0..MMAX] ]);
rkM := Rank(AM);
bE := Matrix(Rationals(), MMAX+1, 1, [ Coefficient(E, m) : m in [0..MMAX] ]);
printf "E in M: %o\n", Rank(HorizontalJoin(AM, bE)) eq rkM;


pairs := [];
for Dp in [ d : d in Divisors(DN) | d gt 1 and IsSquarefree(d) and IsOdd(#PrimeDivisors(d)) ] do
    for Rl in Divisors(DN div Dp) do
        Append(~pairs, <Dp, Rl>);
    end for;
end for;
printf "definite structures: %o\n", pairs;
dS := Dimension(S);
printf "GUARD: DEPTH = %o must exceed ngross + dim S = %o + %o = %o\n",
    DEPTH, #pairs, dS, #pairs + dS;
error if DEPTH le #pairs + dS,
    "DEPTH too small: the rank test would be vacuous (more columns than rows)", DEPTH, #pairs + dS;

grosscol := AssociativeArray(); grossname := AssociativeArray();
for pr in pairs do
    Bq := QuaternionAlgebra(pr[1]);
    Oq := MaximalOrder(Bq);
    OR := pr[2] eq 1 select Oq else Order(Oq, pr[2]);
    CM := Matrix(Rationals(), 4, 4, &cat[ Eltseq(x) : x in Basis(OR) ]);
    LZ := Lattice(CM);
    Bvecs := [ Bq! Eltseq(v) : v in Basis(LZ) ];
    den0 := Lcm([Denominator(Trace(x)) : x in Bvecs]);
    TrInt := Matrix(Integers(), #Bvecs, 1, [ Integers()!(den0*Trace(x)) : x in Bvecs ]);
    KintZ := KernelMatrix(TrInt);
    S0 := [ &+[ KintZ[i][j]*Bvecs[j] : j in [1..#Bvecs] ] : i in [1..Nrows(KintZ)] ];
    GG := Matrix(Rationals(), 3,3, [ Trace(S0[i]*Conjugate(S0[j])) : j in [1..3], i in [1..3] ]);
    LGr := LatticeWithGram(GG);
    repsG := Representatives(Genus(LGr));
    num := [ Rationals()| 0 : m in [0..MMAX] ]; den := 0;
    for Lr in repsG do
        wt := 1/#AutomorphismGroup(Lr);
        T := ThetaSeries(Lr, 2*MMAX+1);
        for m in [0..MMAX] do num[m+1] +:= wt*Coefficient(T, 2*m); end for;
        den +:= wt;
    end for;
    grosscol[pr] := [ x/den : x in num ];
    grossname[pr] := Sprintf("GROSS(%o,%o) disc %o, %o cls", pr[1], pr[2], Determinant(GG)/2, #repsG);
    printf "  built %o\n", grossname[pr];
end for;

cuspcols := [ [ Coefficient(BS[j], m) : m in [0..MMAX] ] : j in [1..#BS] ];
v := Vector(Rationals(), [ Coefficient(E, m) : m in [0..MMAX] ]);

runtest := procedure(tag, prs)
    cols := [* *]; names := [];
    for pr in prs do Append(~cols, grosscol[pr]); Append(~names, grossname[pr]); end for;
    for j in [1..#cuspcols] do Append(~cols, cuspcols[j]); Append(~names, Sprintf("CUSP_%o", j)); end for;
    nc := #cols;
    A := Matrix(Rationals(), MMAX+1, nc, [ [cols[j][m+1] : j in [1..nc]] : m in [0..MMAX] ]);
    rA := Rank(A); rAug := Rank(HorizontalJoin(A, Matrix(MMAX+1,1,Eltseq(v))));
    printf "RESULT %o: %o columns, %o rows, rank %o vs %o -> E in span: %o\n",
        tag, nc, MMAX+1, rA, rAug, rA eq rAug;
    if rA eq rAug then
        x := Solution(Transpose(A), v);
        printf "   E_eis =\n";
        for j in [1..nc] do
            if x[j] ne 0 then printf "     %8o * %o\n", x[j], names[j]; end if;
        end for;
        printf "   aliasing dim %o\n", Dimension(Kernel(Transpose(A)));
    end if;
end procedure;

printf "\n==== UNRESTRICTED (all definite structures) ====\n";
runtest("ALL", pairs);

printf "\n==== Q4: the five preregistered candidate supports ====\n";
// every candidate must be a DEFINITE structure that actually got built --
// D' needs an ODD number of primes (a definite quaternion discriminant has an
// even number of ramified places counting infinity) and R | DN/D'.  Check it
// explicitly rather than index into grosscol and get "value not set".
cands := [
    <"385", [ <385,1>, <385,3> ]>,
    <"105", [ <105,1>, <105,11> ]>,
    <"11",  [ <11,1>,  <11,105> ]>,
    <"231", [ <231,1>, <231,5> ]>,
    <"165", [ <165,1>, <165,7> ]>
];
for c in cands do
    for pr in c[2] do
        error if not IsDefined(grosscol, pr),
            "candidate support is not a built definite structure", c[1], pr;
    end for;
end for;
for c in cands do runtest(c[1], c[2]); end for;
printf "DONE\n";
quit;
