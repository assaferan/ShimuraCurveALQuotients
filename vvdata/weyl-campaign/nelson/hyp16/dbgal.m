AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/nelson/nelson3.m";
D := 2002; N := 1; PREC := 30; CC := ComplexField(PREC); RR := RealField(PREC); Ic := CC.1;
Wtab := WkTable(2, 22, RR);
dat := NelsonSetupT(D, N, 2, PREC, 22, Wtab);
A := dat`A; O := dat`O; mrep := dat`mrep;
bas := [ A ! b : b in ZBasis(O) ];
printf "basis norms %o traces %o\n", [ Norm(b) : b in bas ], [ Trace(b) : b in bas ];
Mb := [ ChangeRing(mrep(b), RR) : b in bas ];
GF4 := Matrix(RR, 4, 4, [ &+[ Mb[i][r][c]*Mb[j][r][c] : r,c in [1,2] ] : i,j in [1..4] ]);
printf "Frobenius Gram determinant %o, eigen-ish diag %o\n", Determinant(GF4), [ GF4[i][i] : i in [1..4] ];
SCL := 10^12;
GI := Matrix(Integers(), 4, 4, [ Round(SCL*GF4[i][j]) : i,j in [1..4] ]);
printf "integral Gram positive definite: %o\n", IsPositiveDefinite(GI);
L := LatticeWithGram(GI);
for m in [11, 91] do
  for bnd in [40, 100, 400] do
    sv := ShortVectors(L, Round(SCL*bnd));
    cnt := 0; nrm := 0; hits := [];
    for t in sv do
      v := Eltseq(t[1]); x := &+[ v[i]*bas[i] : i in [1..4] ];
      if Trace(x) ne 0 then continue; end if;
      cnt +:= 1;
      if Norm(x) ne m then continue; end if;
      nrm +:= 1;
      M := ChangeRing(mrep(x), RR);
      a := M[1][1]; c := M[2][1];
      if Abs(c) lt 1e-12 then Append(~hits, <"c=0", 0>); continue; end if;
      z := (a + Sqrt(RR!m)*Ic)/c; if Imaginary(z) lt 0 then z := (a - Sqrt(RR!m)*Ic)/c; end if;
      Append(~hits, <ChangeRing(z, ComplexField(8)), Round(1000*Imaginary(z))/1000>);
    end for;
    printf "m=%o bnd=%o: %o short vectors, %o with trace 0, %o with norm %o; first fixed points %o\n",
       m, bnd, #sv, cnt, nrm, m, hits[1..Minimum(6,#hits)];
  end for;
end for;
printf "DBGAL DONE\n"; quit;
