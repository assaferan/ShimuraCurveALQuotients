// Robustness check: is #L_quo (the candidate "glue index") stable across order representatives?
// If it wobbles, it is not a well-defined invariant and cannot be the m=0 glue factor.
AttachSpec("ShimuraQuotients.spec");
fh := Open("repcheck_out.txt", "w");
emit := procedure(s) fprintf fh, "%o\n", s; Flush(fh); end procedure;

glue_index := function(Ld, d)
    Q := ChangeRing(Ld`Q, Integers());
    lam := ElementsOfNorm(Ld`Q, [-d], Ld`O, Ld`basis_L)[-d];
    c := Content(lam);
    Lplus := RSpaceWithBasis(Matrix(lam div c));
    Lminus := Kernel(Transpose(Matrix(lam*Q)));
    L := RSpaceWithBasis(IdentityMatrix(Integers(),3));
    return #(L / (Lplus + Lminus));
end function;

D, N := 15, 2;
Ld1 := ShimuraCurveLattice(D, N);
// find a genuinely different order representative (different Gram), as SchoferIsometry.m does
found := false; Ld2 := Ld1;
for tries in [1..80] do
    _ := [Random(1, 10^6) : i in [1..(tries mod 11) + 1]];
    cand := ShimuraCurveLattice(D, N);
    if Eltseq(cand`Q) ne Eltseq(Ld1`Q) then found := true; Ld2 := cand; break; end if;
end for;
emit(Sprintf("found different rep: %o", found));
emit(Sprintf("Q1 = %o", Eltseq(Ld1`Q)));
emit(Sprintf("Q2 = %o", Eltseq(Ld2`Q)));
for d in [-7,-15,-60,-52,-88] do
    g1 := glue_index(Ld1, d);
    g2 := glue_index(Ld2, d);
    emit(Sprintf("  d=%-5o  glue(rep1)=%-4o  glue(rep2)=%-4o  %o", d, g1, g2, g1 eq g2 select "STABLE" else "WOBBLES"));
end for;
emit("DONE");
