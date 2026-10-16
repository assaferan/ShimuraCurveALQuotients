// Is the refutation caused by bad DATA or by a bad DERIVATION?  Same three independent checks
// used at X_0^6(379), now at the two genus-2 star curves where the criterion overshoots worst.
// If the data checks out here too, the criterion's derivation is what is wrong.
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
for b in [ <6,403>, <3390,1>, <2418,1> ] do
    D := b[1]; N := b[2]; DN := D*N;
    full := { Integers() | d : d in Divisors(DN) | GCD(d, DN div d) eq 1 };
    X := CreateShimuraQuot(D, N, full); X`g := GenusShimuraCurveQuotient(D, N, full);
    V := ArithmeticVolume(FuchsianGroup(MaximalOrder(QuaternionAlgebra(D))));
    psi := N; for pr in PrimeDivisors(N) do psi := psi*(pr+1) div pr; end for;
    vq := V*psi/#full;
    ell := NumberOfEllipticPointsByCMOrder(X);
    printf "\n=== X_0^%o(%o)*  genus %o  volume %o\n", D, N, X`g, vq;
    s := vq; tot := 0;
    for q in Sort(SetToSequence(Keys(ell))) do
        for d in Sort(SetToSequence(Keys(ell[q]))) do
            n := ell[q][d]; if n eq 0 then continue; end if;
            tot +:= n; s -:= n*(1-1/q);
            ok := true; deg := 0;
            try deg := Degree(AbsoluteField(FieldOfDefinitionOfCMPoint(X,d))); catch e ok:=false; end try;
            h := ClassNumber(QuadraticForms(d));
            printf "   q=%o d=%-8o n=%-3o h(d)=%-5o deg=%-4o h/deg=%o\n", q, d, n, h,
               ok select Sprint(deg) else "n/a",
               (ok and deg gt 0 and h mod deg eq 0) select
                  Sprintf("%o%o", h div deg, (#full mod (h div deg) eq 0) select " ok" else " <== does not divide #W*")
               else "NOT INTEGRAL <== PROBLEM";
        end for;
    end for;
    printf "   total elliptic points %o; genus from the volume = %o (computed %o): %o\n",
       tot, 1 + s/2, X`g, (1 + s/2) eq X`g select "CONSISTENT" else "INCONSISTENT <== data problem";
end for;
printf "CHKG2 DONE\n";
quit;
