AttachSpec("ShimuraQuotients.spec");
fh := Open("sig_out.txt","w");
emit := procedure(s) fprintf fh, "%o\n", s; Flush(fh); end procedure;
for base in [<15,2>,<21,2>,<10,11>] do
    D,N := Explode(base);
    Ld := ShimuraCurveLattice(D,N);
    Q := ChangeRing(Ld`Q, Rationals());
    // signature via eigenvalue signs of the Gram (Q_code) and of -Q_code
    ev := [x[1] : x in Roots(CharacteristicPolynomial(ChangeRing(Q,RealField(30))))];
    pos := #[e : e in ev | e gt 0]; neg := #[e : e in ev | e lt 0];
    emit(Sprintf("X0^%o(%o): det(Q_code)=%o", D, N, Determinant(ChangeRing(Ld`Q,Integers()))));
    emit(Sprintf("   Q_code  signature (b+,b-) = (%o,%o)   => kappa = 3/2 - %o/2 + %o/2 = %o",
                 pos, neg, neg, pos, 3/2 - neg/2 + pos/2));
    emit(Sprintf("  -Q_code  signature (b+,b-) = (%o,%o)   => kappa = 3/2 - %o/2 + %o/2 = %o",
                 neg, pos, pos, neg, 3/2 - pos/2 + neg/2));
end for;
emit("DONE");
