SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/nelson/nelson.m";

D := 6; N := 1; k := 4; prec := 40;
CC := ComplexField(prec); RR := RealField(prec);

M := CuspForms(Gamma0(6), 4);
NF := Newforms(M);
printf "#newform orbits at level 6 weight 4: %o\n", #NF;
f := NF[1][1];
printf "f = %o\n", f;
nmax := 4000;
q := qExpansion(f, nmax+1);
afs := [ Rationals() ! Coefficient(q, n) : n in [1..nmax] ];
printf "a_f(1..12) = %o\n", afs[1..12];
af := func< n | n le nmax select CC ! afs[n] else 0 >;

dat := NelsonData(D, N, k, prec);
printf "DN = %o, divisors %o\n", dat`DN, dat`ds;
printf "trace gram det = %o\n", Determinant(dat`T);
for d in dat`ds do
  B := Matrix(Rationals(), 4, 4, [ 0 : i in [1..16] ]);
  // index check: [R^(d):R] should be d^2
  bs := dat`lat[d];
  // express bs in terms of dat`bas
  cols := [];
  for b in bs do
    Append(~cols, Eltseq(b));
  end for;
  // coordinates of bas
  Bm := Matrix(Rationals(), 4, 4, [ Eltseq(x) : x in dat`bas ]);
  Cm := Matrix(Rationals(), 4, 4, cols);
  idx := 1/Determinant(Cm * Bm^(-1));
  printf "  d=%o  index [R^(d):R] = %o   (want %o)\n", d, idx, d^2;
end for;

// --- diagonal: must be real and positive
for z0 in [ CC ! (0.17 + 1.3*CC.1), CC ! (-0.4 + 0.9*CC.1) ] do
  for cut in [25, 35, 45] do
    S, nt := NelsonPair(dat, z0, z0, af, RR!cut, CC);
    printf "z=%o cut=%o : S = %o   (nterms %o)\n", z0, cut, S, nt;
  end for;
end for;
quit;
