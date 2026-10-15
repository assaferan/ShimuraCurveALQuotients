// Does the Hecke coset search find ALL p+1 cosets?
//
// Gamma \ {x in O : Nrd(x) = p} has exactly p+1 elements for p not dividing D*N, because the
// left ideals of reduced norm p are the p+1 lines of O/pO = M_2(F_p) and every one of them is
// principal (indefinite algebra, one-sided class number 1).  So the count is known in advance and
// a search that returns fewer has missed cosets -- which is what the coordinate box did.
//
// Two things are checked here, not one.  COMPLETENESS: the search returns p+1.  And NO
// OVER-COUNTING: enumerating the whole Frobenius ball at the bound that was needed, and
// classifying every element of reduced norm p in it, still gives exactly p+1 classes, so the
// test used to decide that two elements lie in the same coset is neither too strict nor too
// loose.  The old box search is run alongside at the same bases for comparison.
//
//   magma -b CASES:=22_1:13,26_1:5,6_1:5,6_1:13,15_1:2 PR:=60 heckecos.m < /dev/null
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vw3.m";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vwhecke.m";

PREC := StringToInteger(PR);
BOX := 5;

// the search this replaces: a cube in the coordinates of a basis of O
function BoxReps(O, p, bound)
    bas := Basis(O); A := Algebra(O);
    reps := [];
    for c in CartesianPower([-bound..bound], 4) do
        x := &+[ c[i]*bas[i] : i in [1..4] ];
        if Norm(x) ne p then continue; end if;
        if exists{ y : y in reps | (A ! y) * (A ! x)^(-1) in O } then continue; end if;
        Append(~reps, x);
        if #reps eq p+1 then break; end if;
    end for;
    return reps;
end function;

printf "%-9o %-4o %-6o %-9o %-7o %-9o %o\n", "base", "p", "p+1", "shortvec", "bound", "boxof5", "ball census";
nfail := 0;
for cs in Split(CASES, ",") do
    ab := Split(cs, ":");
    xs := Split(ab[1], "_");
    D := StringToInteger(xs[1]); N := StringToInteger(xs[2]); p := StringToInteger(ab[2]);
    error if (D*N) mod p eq 0, Sprintf("p = %o divides D*N = %o", p, D*N);
    A := QuaternionAlgebra(D); O := QuaternionOrder(A, N);
    mrep := FuchsianMatrixRepresentation(A : Precision := PREC);
    t0 := Cputime();
    reps, c := HeckeReps(O, p : mrep := mrep);
    tsv := Cputime(t0);

    // every representative must really have reduced norm p and lie in O
    for x in reps do
        error if Norm(x) ne p, "a representative has the wrong reduced norm";
        error if not ((A ! x) in O), "a representative is not in O";
    end for;

    // census of the whole ball at the bound that was needed
    bas := Basis(O);
    L := LatticeWithGram(MajorantGram(O, mrep));
    cls := []; ntot := 0;
    for v in ShortVectors(L, SCL*c*p) do
        x := &+[ (Integers() ! v[1][i]) * bas[i] : i in [1..4] ];
        if Norm(x) ne p then continue; end if;
        ntot +:= 1;
        if exists{ y : y in cls | (A ! y) * (A ! x)^(-1) in O } then continue; end if;
        Append(~cls, x);
    end for;

    box := BoxReps(O, p, BOX);
    ok := (#reps eq p+1) and (#cls eq p+1);
    if not ok then nfail +:= 1; end if;
    printf "%-9o %-4o %-6o %-9o %-7o %-9o %o classes from %o elements  %o (%o s)\n",
        ab[1], p, p+1, #reps, c, #box, #cls, ntot, ok select "" else "  *** WRONG ***",
        RealField(3)!tsv;
end for;
printf "HECKECOS DONE: %o failures\n", nfail;
quit;
