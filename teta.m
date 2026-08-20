// Is Magma's DedekindEta accurate to the working precision, including at small Im?
// Reference: eta(z) = q^{1/24} prod (1-q^n),  q = e(z)  -- accurate when Im(z) is not small;
// for small Im use eta(-1/z) = sqrt(z/i) eta(z) to move to a large-Im point first.
prod_eta := function(z)
    CC := Parent(z);
    q := Exp(2*Pi(CC)*CC.1*z);
    // number of terms so that |q|^n < 10^-(prec+10)
    nt := Ceiling((Precision(CC)+15)*Log(10)/(2*Pi(CC)*Im(z))) + 5;
    p := CC!1;
    for k in [1..nt] do p *:= (1 - q^k); end for;
    return Exp(2*Pi(CC)*CC.1*z/24)*p;
end function;

for prec in [30, 60, 120, 230] do
    CC := ComplexField(prec); ii := CC.1;
    printf "--- precision %o ---\n", prec;
    for z in [CC | 0.3 + 2.0*ii, 0.1 + 0.5*ii, 0.37 + 0.05*ii, 0.123 + 0.000277*ii] do
        d := DedekindEta(z);
        if Im(z) ge 0.4 then
            r := prod_eta(z);
        else
            // move by S: eta(z) = eta(-1/w) with w = -1/z, and eta(-1/w) = sqrt(w/i) eta(w)
            w := -1/z;
            r := Sqrt(w/ii)*prod_eta(w);
        end if;
        printf "  z = %-28o |eta| ~ 10^%-8o rel.err = %o\n", ChangePrecision(z,8),
               Round(Log(10, Abs(d))), ChangePrecision(Abs(d-r)/Abs(r), 6);
    end for;
end for;
quit;
