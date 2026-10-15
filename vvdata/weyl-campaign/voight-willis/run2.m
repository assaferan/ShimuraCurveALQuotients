load "/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/vw6/vw2.m";
SetColumns(0);
PREC := StringToInteger(PR); NT := StringToInteger(NTERMS); Q := StringToInteger(QQ);
CC := ComplexField(PREC);
A := QuaternionAlgebra(StringToInteger(DISC)); O := QuaternionOrder(A, 1);
G := FuchsianGroup(O); sig := Signature(G); g0 := sig[1]; es := sig[2];
DD := UnitDisc( : Precision := PREC);
FD := FundamentalDomain(G, DD);
rho := Max([Abs(ComplexValue(v)) : v in FD]);
H := UpperHalfPlane();
p := ComplexValue(DiscToPlane(H, DD ! 0) : Precision := PREC);
mrep := FuchsianMatrixRepresentation(A : Precision := PREC);
printf "D=%o sig (%o; %o) rho=%o prec=%o N=%o Q=%o rho^N=%o\n",
   DISC, g0, es, RealField(6)!rho, PREC, NT, Q, RealField(3)!(rho^NT);
t := Cputime();
samp := ReducedSamples(G, O, DD, mrep, p, rho, Q, CC);
printf "samples outside the domain: %o of %o  (%o s)\n", #samp, Q, RealField(4)!Cputime(t);
mx := Max([Abs(s[2]) : s in samp]);
printf "max |w'| over samples: %o\n", RealField(5)!mx;
printf "|jD| over samples: min %o max %o\n", RealField(4)!Min([Abs(DiscJ(s[3],s[1])) : s in samp]), RealField(4)!Max([Abs(DiscJ(s[3],s[1])) : s in samp]);
function DimMk(k) return (k-1)*(g0-1) + &+[ Floor((k/2)*(1 - 1/e)) : e in es ]; end function;
nok := 0; ntot := 0;
for k in [2, 4, 6, 8, 10, 12, 14] do
    exp := DimMk(k); if k eq 2 then exp := Max(exp, g0); end if;
    ntot +:= 1;
    sv := SingVals(NormaliseColumns(RowNormalise(AutSystemNorm(samp, k, NT, p, CC))));
    best := 0.0; bi := 0;
    for i in [1..Min(10,#sv)-1] do
        if sv[i] gt 0 and sv[i+1]/sv[i] gt best then best := sv[i+1]/sv[i]; bi := i; end if;
    end for;
    ok := (bi eq exp); if ok then nok +:= 1; end if;
    printf "k=%2o expect %o got %o %-8o gapratio %o  sv[1..5] %o\n", k, exp, bi,
        ok select "OK" else "MISMATCH", RealField(4)!best, [RealField(3)!s : s in sv[1..Min(5,#sv)]];
end for;
printf "RUN2: %o of %o weights correct\n", nok, ntot;
quit;
