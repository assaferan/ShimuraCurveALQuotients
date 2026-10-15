// Independent check of every equation in equations.m: the genus from the genus formula, the genus
// of the curve itself, and the numbers of points over F_p and F_{p^2} against the trace formula at
// every odd prime of good reduction below 60.  Nothing here uses any step of the construction.
//   magma -b logfile:=... verify.m < /dev/null

AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
L := Open(logfile, "w");
procedure Say(s) fprintf L, "%o\n", s; Flush(L); end procedure;

load "vvdata/weyl-campaign/molina/equations.m";

ntot := 0;
for row in molina_reproduced cat new_equations do
    D := row[1]; m := row[2]; g := row[3]; R := row[4];
    W := (m eq 1) select {Integers() | 1} else {Integers() | 1, m};
    X := CreateShimuraQuot(D, 1, W);
    X`g := g;
    okg := GenusShimuraCurveQuotient(D, 1, W) eq g;
    C := HyperellipticCurve(R);
    okc := Genus(C) eq g;
    nok := 0; bad := [];
    for p in PrimesUpTo(60) do
        if p eq 2 or D mod p eq 0 then continue; end if;
        fp := ChangeRing(R, GF(p));
        if Discriminant(fp) eq 0 or Degree(fp) lt 2*g+1 then continue; end if;
        for d in [1, 2] do
            a := #Points(ChangeRing(C, GF(p^d)));
            b := ComputePointsViaTrace(X, p, d);
            if a eq b then nok +:= 1; else Append(~bad, <p, d, a, b>); end if;
        end for;
    end for;
    ntot +:= nok;
    Say(Sprintf("D = %o, w_%o: genus formula %o, curve genus %o, %o point counts agree, %o disagree%o",
        D, m, okg select "ok" else "WRONG", okc select "ok" else "WRONG", nok, #bad,
        #bad eq 0 select "" else Sprintf(" %o", bad)));
end for;
Say(Sprintf("TOTAL %o point counts against the trace formula", ntot));
Say("VERIFY DONE");
quit;
