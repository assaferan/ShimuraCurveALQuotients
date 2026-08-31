// Does IntegralSolution rescue the NONINTEGRAL class?
//
//   magma -b DD:=33 NN:=2 intsol.m
//
// Runs BorcherdsForms twice on the same base -- default, then IntegralSolution -- and reports
// for each whether it succeeded and whether the resulting infinity-side principal parts are
// INTEGRAL.  The point is the pair: a base is "rescued" iff default gives NONINTEGRAL forms (or
// a wrong/inconsistent answer) and IntegralSolution gives integral ones.
//
// Non-integrality is read off the same way ppint.m does it: the principal part of each form at
// infinity must have integral coefficients (Schofer's c_eta(m) in Z).
AttachSpec("ShimuraQuotients.spec");
D := 33; N := 2;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;
printf "INTSOL BASE %o_%o\n", D, N;

Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};

// principal-part integrality of one form at infinity
function ppIntegral(f)
    q := qExpansionAtoo(f, 1);
    v := Valuation(q);
    if v ge 0 then return true, 1; end if;
    dens := [Denominator(Rationals()!Coefficient(q, i)) : i in [v..-1]];
    if IsEmpty(dens) then return true, 1; end if;
    L := LCM(dens);
    return L eq 1, L;
end function;

procedure run(star, curves, D, N, useint)
    tag := useint select "INTSOL" else "DEFAULT";
    t0 := Realtime();
    try
        fs := BorcherdsForms(star, curves : Prec := 100, IntegralSolution := useint);
        allint := true; maxden := 1;
        for k in Keys(fs) do
            ok, L := ppIntegral(fs[k]);
            if not ok then allint := false; end if;
            if L gt maxden then maxden := L; end if;
        end for;
        printf "INTSOL %o_%o %o OK nforms %o integral %o maxden %o time %os\n",
               D, N, tag, #Keys(fs), allint, maxden, Realtime() - t0;
    catch e
        obj := Sprint(e`Object);
        short := (#obj gt 90) select obj[1..90] else obj;
        printf "INTSOL %o_%o %o FAIL time %os :: %o\n", D, N, tag, Realtime() - t0,
               Join(Split(short, "\n"), " | ");
    end try;
end procedure;

run(star, curves, D, N, false);
run(star, curves, D, N, true);
printf "INTSOL DONE %o_%o\n", D, N;
quit;
