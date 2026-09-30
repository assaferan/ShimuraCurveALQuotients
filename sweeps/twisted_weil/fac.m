// Factor the hit polynomials.  From the repo root: magma -b [IN:=...] sweeps/twisted_weil/fac.m < /dev/null
SetColumns(0); if not assigned IN then IN := "sweeps/twisted_weil/hitpolys.txt"; end if; _<t> := PolynomialRing(Integers());
for l in Split(Read(IN), "\n") do f := Split(l, " "); if #f lt 2 then continue; end if;
  c := eval f[2]; P := &+[c[i]*t^(#c-i) : i in [1..#c]];
  printf "p=%o %o = %o\n", f[1], f[2], Factorization(P);
end for; quit;
