// deg Z(d) -- the DEGREE of the disc-d CM cycle on the star curve, i.e. the number of points
// counted over Qbar.  This is the quantity that prices an even correction: perturbing the divisor
// by <d, amt> raises deg f by amt * deg Z(d), and RationalConstraintsOnEquations then needs
// 2g+5 + amt*deg Z(d) RATIONAL CM points to fit it (QUADCONSTRAINTS.md sec 8).
//
// deg Z(d) = sum of the degrees of the fields of definition -- the same quantity replace_column
// (SchoferFormula.m:1740) already uses as `deg`.
//
// VALIDATION IS BUILT IN: at 34_3 the divisor-degree identity (QUADCONSTRAINTS.md sec 7b) pins
// deg Z(3) = deg Z(24) = deg Z(51) = deg Z(408) = 1 and deg Z(68) = 2, independently of anything
// here.  Those five are checked before any new value is printed.
//
//   magma -b D_s:=34 N_s:=3 vvdata/weyl-campaign/even-correction/degz.m < /dev/null
SetQuitOnError(true);
SetColumns(0);
AttachSpec("ShimuraQuotients.spec");
D := StringToInteger(D_s); N := StringToInteger(N_s);
bd := 600; if assigned BD then bd := StringToInteger(BD); end if;
curves := GetHyperellipticCandidates();
Xstar := rep{X : X in curves | X`D eq D and X`N eq N and IsStarCurve(X)};

function degZ(X, d)
    try
        flds := FieldsOfDefinitionOfCMPointFast(X, d);
    catch e
        return -1;
    end try;
    if IsEmpty(flds) then return 0; end if;
    return &+[Degree(F) : F in flds];
end function;

// ---- validation against the divisor-degree identity, before trusting anything new ----
known := [<3,1>, <24,1>, <51,1>, <408,1>, <68,2>];
nok := 0;
printf "VALIDATION (expected from the divisor-degree identity at 34_3):\n";
for t in known do
    got := degZ(Xstar, -t[1]);
    ok := (got eq t[2]);
    if ok then nok +:= 1; end if;
    printf "  deg Z(%o) = %o   expected %o   %o\n", -t[1], got, t[2], ok select "OK" else "*** MISMATCH ***";
end for;
printf "  => %o of %o known values reproduced\n\n", nok, #known;
if (D eq 34) and (N eq 3) and (nok ne #known) then
    printf "REFUSING to print a sweep: the method does not reproduce the known values.\n";
    quit;
end if;

printf "SWEEP deg Z(d) for D=%o N=%o, |d| <= %o\n", D, N, bd;
for a in [1..bd] do
    d := -a;
    if not (d mod 4 in {0,1}) then continue; end if;
    z := degZ(Xstar, d);
    if z le 0 then continue; end if;
    printf "DEGZ d %o degZ %o\n", d, z;
end for;
printf "DEGZ_DONE\n";
quit;
