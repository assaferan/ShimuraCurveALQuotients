// Gate for etared.m: etaRed at MODEST precision must reproduce DedekindEta
// computed at HUGE precision, on exactly the arguments the 1155 E-table uses.
//   magma -b vvdata/weyl-campaign/etared_check.m
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/fastcosets.m";
load "vvdata/weyl-campaign/etared.m";

M := 4620; ds := Divisors(M);
reps := fastCosetReps(M); words := [ VVSTWord(g) : g in reps ];

REF := ComplexField(2000);        // ground truth
tau0R := REF!0.31 + REF!1.31*REF.1;
slash := function(word, tau, CC)
    z := tau; factor := CC!1;
    for i := #word to 1 by -1 do
        if word[i][1] eq "S" then factor /:= Sqrt(z); z := -1/z; else z := z + word[i][2]; end if;
    end for;
    return factor, z;
end function;

for PREC in [120, 200] do
    CC := ComplexField(PREC);
    tau0 := CC!0.31 + CC!1.31*CC.1;
    worst := RealField(20)!0; wat := "";
    // (a) the sfun arguments tau0/d at the S-coset: (a,b,e) = (1,0,d)
    for d in ds do
        v := etaRed(tau0/d);
        ref := DedekindEta(tau0R/d);
        rel := Abs(REF!v - ref)/Abs(ref);
        if rel gt worst then worst := RealField(20)!rel; wat := Sprintf("tau0/%o", d); end if;
    end for;
    // (b) the num arguments d*z0 at a spread of cosets, including the costly
    //     c = 1 block where Im(z0) ~ 1e-7
    for wi in [3, 500, 1000, 2000, 2500, 5000, 10000, 13800] do
        _, z0 := slash(words[wi], tau0, CC);
        _, z0R := slash(words[wi], tau0R, REF);
        for d in [1, 2, 7, 105, 660, 4620] do
            v := etaRed(d*z0);
            ref := DedekindEta(d*z0R);
            rel := Abs(REF!v - ref)/Abs(ref);
            if rel gt worst then worst := RealField(20)!rel; wat := Sprintf("coset %o, d = %o", wi, d); end if;
        end for;
    end for;
    printf "PREC %o: worst RELATIVE |etaRed - DedekindEta@2000| = %o   (at %o)\n",
        PREC, worst, wat;
    printf "PREC %o VERDICT: %o\n", PREC,
        worst lt (RealField(20)!10)^(-(PREC - 20)) select "PASS" else "FAIL";
end for;

// speed: the whole point
CC := ComplexField(120);
tau0 := CC!0.31 + CC!1.31*CC.1;
for wi in [3, 1000, 2500] do
    _, z0 := slash(words[wi], tau0, CC);
    t := Cputime(); for r in [1..20] do x := etaRed(z0); end for; tr := Cputime(t)/20;
    t := Cputime(); for r in [1..3] do x := DedekindEta(z0); end for; td := Cputime(t)/3;
    printf "coset %5o: etaRed %o ms   DedekindEta %o ms\n", wi,
        RealField(6)!(1000*tr), RealField(6)!(1000*td);
end for;
printf "DONE\n";
quit;
