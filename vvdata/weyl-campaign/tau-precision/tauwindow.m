// REPAIR track, closing half: WHY does wi = 2 fail, with a healthy Im(z0) = 0.72?
//
// eta is never evaluated at z0 itself.  In M0MultiplierExact:
//    num0        = fac0 * prod_i DedekindEta(d_i * z0)^r_i           d_i in ds = Divisors(M)
//    sfun(tau0)  = ee(...) * prod_i DedekindEta((a_i*tau0+b_i)/e_i)^r_i   from triang(g, d)
//    k0          = num0 / sfun(tau0)
//
// For Im(z) large, |eta(z)| ~ |q^(1/24)| = exp(-2*pi*Im(z)/24), so
//    log10|eta(z)| ~ -(2*pi/(24*ln10)) * Im(z) = -0.11374 * Im(z).
// With d up to M = 1220, a LARGE Im(z0) makes Im(d*z0) enormous and eta underflow to nothing;
// a TINY Im(z0) does the opposite.  Both extremes destroy the ratio in floating point even
// though the ratio is mathematically fine.  Predicted: there is a SAFE WINDOW of Im(z0), and
// wi = 2 sits above it while wi = 1221 sits below it.
AttachSpec("ShimuraQuotients.spec");

D := 10; N := 61;
M := IsOdd(D*N) select 4*D*N else 2*D*N;
ds := Divisors(M);
printf "M = %o   ds = %o\n\n", M, ds;

CC := ComplexField(80); ii := CC.1;
tau0 := CC!0.31 + CC!1.31*ii;

slashdata := function(word, tau)
    z := tau; factor := CC!1;
    for i := #word to 1 by -1 do
        if word[i][1] eq "S" then factor /:= Sqrt(z); z := -1/z;
        else z := z + word[i][2]; end if;
    end for;
    return factor, z;
end function;

// verbatim from VectorValuedForm.m:395-408
triang := function(g, d)
    g2 := Matrix(Integers(), 2, 2, [d*g[1][1], d*g[1][2], g[2][1], g[2][2]]);
    c1 := g2[1][1]; c2 := g2[2][1];
    h := GCD(c1, c2);
    p1 := c1 div h; p2 := c2 div h;
    gg, u, v := XGCD(p1, p2);
    error if gg ne 1, "triangularisation: row not primitive";
    gd := Matrix(Integers(), 2, 2, [p1, -v, p2, u]);
    sd := gd^(-1) * g2;
    a := sd[1][1]; b := sd[1][2]; e := sd[2][2];
    if a lt 0 then a := -a; b := -b; e := -e; end if;
    if e lt 0 then gd := -gd; sd := -sd; a := sd[1][1]; b := sd[1][2]; e := sd[2][2]; end if;
    return a, b, e;
end function;

reps := VVCosetReps(M);
words := [ VVSTWord(g) : g in reps ];
printf "#words = %o\n\n", #words;

C := 2*Pi(RealField(30))/(24*Log(RealField(30)!10));   // 0.11374...

// log10|eta| at height h, via whichever cusp expansion converges:
//   large Im: -C*Im.   small Im: use eta(-1/z) = sqrt(-iz) eta(z) => |eta(z)| ~ |z|^(-1/2)*...
function logeta(z)
    h := Imaginary(z);
    if h ge 0.05 then return -C*h; end if;
    w := -1/z;                                   // Im(w) = Im(z)/|z|^2, large
    return -0.5*Log(Abs(z))/Log(10.0) - C*Imaginary(w);
end function;

procedure report(wi, words, ds, tau0, C)
    _, z0 := slashdata(words[wi], tau0);
    g := VVWordMatrix(words[wi]);
    lo_n := 10^30; hi_n := -10^30; lo_s := 10^30; hi_s := -10^30;
    for d in ds do
        ln := logeta(d*z0);
        if ln lt lo_n then lo_n := ln; end if;
        if ln gt hi_n then hi_n := ln; end if;
        a, b, e := triang(g, d);
        ls := logeta((a*tau0 + b)/e);
        if ls lt lo_s then lo_s := ls; end if;
        if ls gt hi_s then hi_s := ls; end if;
    end for;
    RR := RealField(6);
    printf "  wi %-6o Im(z0) %-12o | log10|eta(d*z0)| in [%o, %o] | triang in [%o, %o] | SPREAD %o digits\n",
        wi, RR!Imaginary(z0), RR!lo_n, RR!hi_n, RR!lo_s, RR!hi_s,
        RR!Maximum(hi_n - lo_n, hi_s - lo_s);
end procedure;

printf "THE TWO FAILING COSETS:\n";
for wi in [2, 1221] do report(wi, words, ds, tau0, C); end for;

printf "\nTYPICAL COSETS (should be far tamer):\n";
for wi in [37, 300, 611, 900, 1500, 2000] do
    if wi le #words then report(wi, words, ds, tau0, C); end if;
end for;

// Where is the safe window?  Scan the whole range and report the worst spread per Im(z0) band.
printf "\nSPREAD vs Im(z0), over all %o cosets:\n", #words;
bands := [ <0.0, 1e-5>, <1e-5, 1e-4>, <1e-4, 1e-3>, <1e-3, 1e-2>,
           <1e-2, 1e-1>, <1e-1, 1.0>, <1.0, 100.0> ];
cnt := [ 0 : b in bands ]; wor := [ RealField(30)!0 : b in bands ];
for wi in [1..#words] do
    _, z0 := slashdata(words[wi], tau0);
    g := VVWordMatrix(words[wi]);
    lo := 10^30; hi := -10^30;
    for d in ds do
        ln := logeta(d*z0);
        if ln lt lo then lo := ln; end if;
        if ln gt hi then hi := ln; end if;
        a, b, e := triang(g, d);
        ls := logeta((a*tau0 + b)/e);
        if ls lt lo then lo := ls; end if;
        if ls gt hi then hi := ls; end if;
    end for;
    h := Imaginary(z0);
    for j->bd in bands do
        if h ge bd[1] and h lt bd[2] then
            cnt[j] +:= 1;
            if hi - lo gt wor[j] then wor[j] := hi - lo; end if;
            break;
        end if;
    end for;
end for;
RR := RealField(6);
printf "   %-22o %-8o %o\n", "Im(z0) band", "#cosets", "worst spread (digits)";
for j->bd in bands do
    printf "   [%-8o, %-8o) %-8o %o\n", RR!bd[1], RR!bd[2], cnt[j], RR!wor[j];
end for;
quit;
