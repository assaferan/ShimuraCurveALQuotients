// What does a base cost once it resolves?  Only two numbers decide it: the time to build the
// fundamental domain, and its circumradius rho, which fixes the term count through
// NT = T ln10 / (-ln rho).  The singular value decomposition is then the whole expense, and it
// scales as NT^3 at fixed working precision -- the precision does not have to grow with NT,
// because the dynamic range of the matrix is rho^NT, which the target T pins to 1e-T by
// construction.  Measured at D = 15, 22 and 26 the constant is 3.1e-8 seconds at precision 200
// for a system of 2*NT rows, and memory for the real form of the matrix is 4*2*NT^2 complex
// entries, about 332 bytes each at that precision.
//
// rho grows towards 1 with the arithmetic volume: a Dirichlet domain of hyperbolic area 2 pi V
// cannot fit inside |w| < r unless r^2/(1-r^2) >= V/2, so rho >= sqrt(V/(V+2)), and the measured
// values run about 1.1 to 1.25 times that.  Since the term count goes as 1/(-ln rho), this is
// what sets the reach.
//
//   magma -b BASES:=330_1,510_1,462_1,6_35 PR:=120 RC:=3 TGT:=20 rhoprobe.m < /dev/null
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vw3.m";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vwal.m";

PREC := StringToInteger(PR); rcit := StringToInteger(RC);
TGT := StringToInteger(TGT);
R4 := RealField(4); RR := RealField(PREC);
H := UpperHalfPlane();

printf "base | volume | domain s | sides | rho | NT at 1e-%o | svd estimate | memory estimate\n", TGT;
for b in Split(BASES, ",") do
    xs := Split(b, "_");
    D := StringToInteger(xs[1]); N := StringToInteger(xs[2]);
    A := QuaternionAlgebra(D); O := QuaternionOrder(A, N);
    G := FuchsianGroup(O);
    t0 := Cputime();
    if rcit gt 0 then _, DD, FD := RecentredDisc(G, H, PREC, rcit);
    else DD := UnitDisc( : Precision := PREC); FD := FundamentalDomain(G, DD); end if;
    td := Cputime(t0);
    rho := Max([Abs(ComplexValue(v)) : v in FD]);
    NT := TermsFor(rho, TGT);
    svd := 3.1e-8 * NT^3;                       // seconds, calibrated at D = 15, 22, 26
    mem := 4.0 * 2 * NT^2 * 332 / 2^30;         // gigabytes for the real form of the matrix
    printf "%o | %o | %o | %o | %o | %o | %o | %o GB\n",
        b, R4!ArithmeticVolume(G), R4!td, #FD, R4!rho, NT,
        svd lt 7200 select Sprintf("%o s", Round(svd)) else Sprintf("%o h", R4!(svd/3600)),
        R4!mem;
end for;
printf "RHOPROBE DONE\n";
quit;
