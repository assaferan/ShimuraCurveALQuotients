SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/nelson/nelson.m";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/nelson/wk.m";
prec := 40; RR := RealField(prec);
for k in [2,4] do
  Wk := func< x | WkFun(k, x, RR) >;
  printf "k=%o: W_k at 1,3,6,12,25,50,80: %o\n", k, [RR!Wk(RR!x) : x in [1,3,6,12,25,50,80]];
  t := Cputime(); P := WkPanels(Wk, RR!0.5, RR!90, RR : ratio := 1.5, deg := 44);
  printf "  panel build %o s (%o panels)\n", RR!Cputime(t), #P[2];
  worst := RR!0;
  for j in [1..40] do
    x := RR ! (0.5 + 89.5*j/41);
    e := Abs(WkFast(P,x) - Wk(x))/Abs(Wk(x));
    worst := Max(worst, e);
  end for;
  printf "  worst relative error of the approximation: %o\n", worst;
  t := Cputime(); for j in [1..2000] do _ := WkFast(P, RR!(1 + 80*j/2000)); end for;
  tf := Cputime(t);
  t := Cputime(); for j in [1..100] do _ := Wk(RR!(1 + 80*j/100)); end for;
  td := Cputime(t);
  printf "  2000 fast evals %o s  vs  100 direct %o s   => speedup %o x\n", RR!tf, RR!td, RR!((td/100)/(tf/2000));
end for;
quit;
