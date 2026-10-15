// One (D,N) per invocation, so a kill or timeout truncates cleanly and does not lose the others.
SetColumns(0);
D := StringToInteger(DD); N := StringToInteger(NN);
printf "D=%o N=%o DN=%o: ", D, N, D*N;
t := Cputime();
A := QuaternionAlgebra(D); O := QuaternionOrder(A, N); G := FuchsianGroup(O);
printf "group %.1o s, ", Cputime(t);
t := Cputime(); P := FundamentalDomain(G); tf := Cputime(t);
printf "domain %o sides %.1o s, ", #P, tf;
t := Cputime(); Gp := Group(G);
printf "presentation %o gens %.1o s\n", #Generators(Gp), Cputime(t);
quit;
