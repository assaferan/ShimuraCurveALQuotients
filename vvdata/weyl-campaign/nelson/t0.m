SetColumns(0);
prec := 40;
RR := RealField(prec); CC := ComplexField(prec);
// KBessel sanity: K_0(1)=0.421024438240708333..., K_1(1)=0.601907230197234574...
printf "K_0(1) = %o\n", KBessel(0, RR!1);
printf "K_1(1) = %o\n", KBessel(1, RR!1);
printf "K_3(12.6) = %o\n", KBessel(3, RR!12.6);
printf "K_2(0.5) = %o   (expect 7.550183551)\n", KBessel(2, RR!0.5);

A := QuaternionAlgebra(6);
O := MaximalOrder(A);
mrep := FuchsianMatrixRepresentation(A : Precision := prec);
printf "mrep of 1: %o\n", mrep(A!1);
for b in Basis(O) do
  M := mrep(A!b);
  printf "b=%o  Nrd=%o  det(iota b)=%o  trd=%o tr(M)=%o\n", b, Norm(A!b), Determinant(M), Trace(A!b), M[1][1]+M[2][2];
end for;
printf "BaseRing mrep image: %o\n", BaseRing(mrep(A!1));
quit;
