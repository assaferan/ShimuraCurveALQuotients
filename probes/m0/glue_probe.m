// Probe: d-dependent glue index |L/(L_+ + L_-)| for the m=0 term (gap 2).
// L = Z^3 with Gram Q; L_+ = Z*(lambda/content); L_- = lambda^perp cap L.
// Reproduces exactly the L_quo that Kappa builds (SchoferFormula.m:591-595), surfaced as a number.
// Writes incrementally to glue_out.txt (Flush each line) so progress is visible under Magma -b buffering.

AttachSpec("ShimuraQuotients.spec");
fh := Open("glue_out.txt", "w");
emit := procedure(s) fprintf fh, "%o\n", s; Flush(fh); end procedure;

glue_index := function(Ld, d)
    Q := ChangeRing(Ld`Q, Integers());
    lam := ElementsOfNorm(Ld`Q, [-d], Ld`O, Ld`basis_L)[-d];
    c_Lplus := Content(lam);
    Lplus := RSpaceWithBasis(Matrix(lam div c_Lplus));
    Lminus := Kernel(Transpose(Matrix(lam*Q)));
    L := RSpaceWithBasis(IdentityMatrix(Integers(),3));
    L_quo := L / (Lplus + Lminus);
    return #L_quo, Moduli(L_quo), (lam*Q, lam);
end function;

emit("START");
for base in [<15,2>, <21,2>] do
    D, N := Explode(base);
    Ld := ShimuraCurveLattice(D, N);
    emit(Sprintf("=== X0^%o(%o), det Q = %o ===", D, N, Determinant(ChangeRing(Ld`Q, Integers()))));
    ds := (D eq 15) select [-7,-15,-52,-60,-88,-120,-148,-168,-228,-232]
                     else  [-7,-15,-16,-28,-60,-120,-148,-168,-228,-232,-408];
    for d in ds do
        try
            g, mdl, nrm := glue_index(Ld, d);
            emit(Sprintf("  d=%-5o  glue=%-4o  moduli=%o  <lam,lam>=%o (expect %o)", d, g, mdl, nrm, 2*d));
        catch e
            emit(Sprintf("  d=%-5o  ERROR: %o", d, e`Object));
        end try;
    end for;
end for;
emit("DONE");
