SetColumns(0);
for D in [410, 665, 435, 390, 930] do
  printf "D = %-5o factorisation %-22o number of primes %o  %o\n", D,
     Sprint(Factorisation(D)), #PrimeDivisors(D),
     IsEven(#PrimeDivisors(D)) select "valid quaternion discriminant"
       else "NOT a discriminant of an indefinite rational quaternion algebra";
end for;
A := QuaternionAlgebra(410);
printf "\nQuaternionAlgebra(410) has discriminant %o, not 410\n", Discriminant(A);
printf "G410B DONE\n"; quit;
