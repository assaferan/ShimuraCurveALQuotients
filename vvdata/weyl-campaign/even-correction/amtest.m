AttachSpec("ShimuraQuotients.spec");
import "TraceFormula.m" : Hurwitz;
function sig(j, p)
    if j lt 0 then return 0; end if;
    return &+[p^i : i in [0..j]];
end function;
function aE(m, D, N)
    D0 := FundamentalDiscriminant(-4*m);
    ok, c := IsSquare((-4*m) div D0);
    if not ok then return 0, false; end if;
    val := 12 * Hurwitz(4*m);
    for p in PrimeDivisors(D) do
        kp := Valuation(c, p); chi := KroneckerSymbol(D0, p);
        den := (p-1)*(sig(kp,p) - chi*sig(kp-1,p));
        if den eq 0 then return 0, false; end if;
        val *:= (1 - chi)/den;
    end for;
    if N gt 1 then
        for p in PrimeDivisors(N) do
            kp := Valuation(c, p); chi := KroneckerSymbol(D0, p);
            den := (p^2-1)*(sig(kp,p) - chi*sig(kp-1,p));
            if den eq 0 then return 0, false; end if;
            val *:= (p - chi)*p^kp/den;
        end for;
    end if;
    return val, true;
end function;
cases := [ <15,2,[<2,-1>,<10,1>,<30,0>]>, <6,5,[<10,3/2>,<15,1/2>,<30,3/2>]>,
           <10,3,[<3,1/2>,<12,1/2>,<30,3/2>]>, <21,2,[<2,-2>,<6,-4>,<18,2>,<42,-8>]> ];
printf "base     m     A_m(meas)   -a_E(m)      -a_E/4       +a_E/4       verdict\n";
nm := 0; nt := 0;
for cs in cases do
  D := cs[1]; N := cs[2];
  for pr in cs[3] do
    m := pr[1]; meas := pr[2]; nt +:= 1;
    v, ok := aE(m, D, N);
    if not ok then printf "%-8o %-5o %-11o ERR\n", Sprintf("%o_%o",D,N), m, meas; continue; end if;
    c1 := -v/4; c2 := v/4;
    tag := (c1 eq meas) select "MATCH -a_E/4" else ((c2 eq meas) select "MATCH +a_E/4" else "no");
    if tag ne "no" then nm +:= 1; end if;
    printf "%-8o %-5o %-11o %-12o %-12o %-12o %o\n", Sprintf("%o_%o",D,N), m, meas, v, c1, c2, tag;
  end for;
end for;
printf "\nMATCHED %o of %o\n", nm, nt;
quit;
