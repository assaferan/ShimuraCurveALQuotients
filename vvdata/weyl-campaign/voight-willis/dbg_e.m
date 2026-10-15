load "/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/vw6/vw2.m";
SetColumns(0);
PREC := 100; CC := ComplexField(PREC);
A := QuaternionAlgebra(6); O := QuaternionOrder(A, 1);
G := FuchsianGroup(O); DD := UnitDisc( : Precision := PREC);
FD := FundamentalDomain(G, DD);
rho := Max([Abs(ComplexValue(v)) : v in FD]);
H := UpperHalfPlane();
p := ComplexValue(DiscToPlane(H, DD ! 0) : Precision := PREC);
mrep := FuchsianMatrixRepresentation(A : Precision := PREC);
samp  := ReducedSamples(G, O, DD, mrep, p, rho, 300, CC);   // fitting set
samp2 := ReducedSamples(G, O, DD, mrep, p, rho*0.97, 97, CC); // FRESH check set, different radius
printf "fit %o samples, check %o samples, rho %o\n", #samp, #samp2, RealField(6)!rho;
// smallest right singular vector of the complex matrix, via the real 2m x 2n form
function SmallestVec(Mx)
    CC := BaseRing(Mx); RR := RealField(Precision(CC));
    m := Nrows(Mx); n := Ncols(Mx);
    MR := Matrix(RR, 2*m, 2*n, [ RR | 0 : i in [1..4*m*n] ]);
    for i in [1..m] do for j in [1..n] do
        a := Real(Mx[i][j]); b := Imaginary(Mx[i][j]);
        MR[i][j] := a; MR[i][n+j] := -b; MR[m+i][j] := b; MR[m+i][n+j] := a;
    end for; end for;
    U, S, V := SingularValueDecomposition(MR);
    svs := [ Abs(S[i][i]) : i in [1..Min(Nrows(S), Ncols(S))] ];
    _, idx := Min(svs);
    v := [ CC | V[j][idx] + CC.1*V[n+j][idx] : j in [1..n] ];
    return v, svs[idx], Sort(svs);
end function;
NT := 50; k := 4;
Mx := AutSystem(samp, k, NT, CC);
v, sval, svs := SmallestVec(Mx);
printf "k=%o N=%o smallest sv (unnormalised matrix) %o\n", k, NT, RealField(4)!sval;
printf "|b_n| for n = 0,5,10,...,45: %o\n", [RealField(3)!Abs(v[n+1]) : n in [0,5,10,15,20,25,30,35,40,45]];
printf "max |b_n| %o at n=%o ; |b_last| %o\n", RealField(4)!Max([Abs(x) : x in v]),
   Index([Abs(x) : x in v], Max([Abs(x) : x in v]))-1, RealField(3)!Abs(v[NT]);
// residual of this vector on the FRESH set
res := [];
for s in samp2 do
    w := s[1]; wp := s[2]; jk := DiscJ(s[3], w)^k;
    Append(~res, Abs(&+[ v[n+1]*(wp^n - jk*w^n) : n in [0..NT-1] ]));
end for;
printf "fresh-set residual: max %o mean %o\n", RealField(4)!Max(res), RealField(4)!(&+res/#res);
printf "DBGE DONE\n";
quit;
