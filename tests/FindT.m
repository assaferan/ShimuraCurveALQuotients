// find_t's LP: Magma gives up sporadically, and the old code read that as infeasibility.
//
// ⚠ WHAT THIS TEST IS FOR.  `find_t` minimises the pole order of the auxiliary eta-quotient `t`
// over an integer LP whose variables carry an ARBITRARY lower bound.  On Magma 2.29-7 the solver
// returns `success = 25` ("gave up") for some (M, bound) pairs, in neither direction monotonically:
//
//     M = 1572, bound -1000  fails, while -261 -300 -500 -800 -1500 -3000 all succeed with k = 260
//     M =  732, bound -5000  fails, while -1000 -1024 -2000 all succeed with k = 120
//
// The old `assert success eq 0` therefore reported "Assertion failed" at M = 1572 (X_0^6(131)) and
// M = 1644 (X_0^6(137)), which reads as an obstruction and is nothing of the kind -- both are
// feasible, with the witness this test re-derives.
//
// ⚠ AND WHY IT CHECKS FEASIBILITY, NOT JUST A NUMBER.  Asserting `k = 2N-2` alone would only say
// the solver returned the value it returned last time.  The property that MATTERS is that the
// returned `t` satisfies the constraints -- which is checkable independently of the solver, and is
// what would catch a "success" carrying junk.  Both are checked; the counts below make a silently
// skipped case fail rather than pass.
import "BorcherdsForms.m" : find_t;
_ := ClassNumberLU(-4);       // AttachSpec is lazy; `import` compiles its file NOW, so touch one
                              // intrinsic first or references from other packages go unresolved.

// M = 12N for these bases, and the optimum is k = 2N-2 with eta exponents that do not depend on N.
// 876/1068/1284 solved at the DEFAULT bound before the retry loop existed -- they are the control
// that the default path still returns exactly what it used to.  1572/1644 used to abort.
KNOWN := [* <876, 144, "default bound">, <1068, 176, "default bound">, <1284, 212, "default bound">,
            <1572, 260, "used to abort">, <1644, 272, "used to abort"> *];
ETA := [0,1,0,-2,-3,6,0,-1,0,2,3,-6];

nfeas := 0; nk := 0;
for c in KNOWN do
    M := c[1];
    t, lhs, rhs, n_eq, n_ds := find_t(M);
    e := Eltseq(t);

    // the objective: the last variable is the pole order
    assert e[#e] eq c[2];
    nk +:= 1;

    // the eta exponents are the first n_ds entries and are N-independent across this family
    assert e[1..n_ds] eq ETA;

    // feasibility, checked against find_t's OWN constraint blocks -- independent of the solver
    v := lhs * Matrix(Integers(), Ncols(lhs), 1, e);
    assert &and[v[i,1] eq rhs[i,1] : i in [1..n_eq]];                       // equalities
    assert &and[v[i,1] ge rhs[i,1] : i in [n_eq+1..n_eq+n_ds]];             // >= block
    idx := n_eq + n_ds + 1;
    assert v[idx,1] le rhs[idx,1];                                          // pole bound above
    assert v[idx+1,1] ge rhs[idx+1,1];                                      // pole bound below
    nfeas +:= 1;
end for;

// ⚠ COUNT THE CHECKS.  If find_t starts erroring, or a case stops being reached, this must go red
// rather than green-with-nothing-verified.
assert nk eq #KNOWN;
assert nfeas eq #KNOWN;
printf "find_t: %o optima reproduced, %o solutions verified feasible...", nk, nfeas;
