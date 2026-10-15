// Theorem 3.1 is stated for a newform of level exactly d_B N.  At N > 1 the differentials of
// X_0^D(N)/W are the W-fixed part of the D-new space of level DN, which contains, besides the
// newforms of level exactly DN, oldform copies of newforms of level DM for proper M | N.  This
// asks, base by base, how much of the space Theorem 3.1 as written can supply: it counts the
// W-fixed newforms of level exactly DN against the genus of the quotient.
//
// Atkin-Lehner, weight 2 and DN squarefree: the classical eigenvalue of w_p is -a_f(p), and the
// quaternionic one is minus that at p | D and equal to it at p | N, so the quaternionic
// eigenvalue of w_m is  prod_{p | m, p | D} a_f(p) * prod_{p | m, p | N} (-a_f(p)).
//
//   magma -b BASES:=141_2:1,3,94,282;39_11:1,11,39,429 nold.m < /dev/null
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
for spec in Split(BASES, ";") do
    pp := Split(spec, ":");
    bn := Split(pp[1], "_");
    D := StringToInteger(bn[1]); N := StringToInteger(bn[2]);
    W := { StringToInteger(x) : x in Split(pp[2], ",") };
    g := GenusShimuraCurveQuotient(D, N, W);
    t0 := Cputime();
    NFs := Newforms(CuspForms(Gamma0(D*N), 2));
    cnt := 0;
    for orb in NFs do
        for f in orb do
            if Level(Parent(f)) ne D*N then continue; end if;
            qq := qExpansion(f, 2*Max(PrimeDivisors(D*N)) + 4);
            inv := true;
            for m in W do
                if m eq 1 then continue; end if;
                e := 1;
                for p in PrimeDivisors(m) do
                    e *:= (D mod p eq 0) select Coefficient(qq, p) else -Coefficient(qq, p);
                end for;
                if e ne 1 then inv := false; end if;
            end for;
            if inv then cnt +:= 1; end if;
        end for;
    end for;
    printf "%o_%o  W=%-30o genus %-3o  W-fixed newforms of level exactly %o: %-3o  %o   (%o s)\n",
        D, N, Sort(SetToSequence(W)), g, D*N, cnt,
        cnt eq g select "COMPLETE -- Theorem 3.1 supplies the whole space"
                 else "SHORT by " cat IntegerToString(g - cnt) cat " -- oldforms needed",
        RealField(4)!Cputime(t0);
end for;
printf "NOLD DONE\n";
quit;
