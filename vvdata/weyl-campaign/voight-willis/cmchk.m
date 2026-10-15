load "/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/vw6/vw3.m";
SetColumns(0);
PREC := 140; NT := 80; CC := ComplexField(PREC);
A := QuaternionAlgebra(6); O := QuaternionOrder(A, 1);
G := FuchsianGroup(O); DD := UnitDisc( : Precision := PREC);
FD := FundamentalDomain(G, DD);
rho := Max([Abs(ComplexValue(v)) : v in FD]);
H := UpperHalfPlane();
p := ComplexValue(DiscToPlane(H, DD ! 0) : Precision := PREC);
mrep := FuchsianMatrixRepresentation(A : Precision := PREC);
samp := ReducedSamples(G, O, DD, mrep, p, rho, 300, CC);
F := func< bs, w | &+[ bs[n+1]*w^n : n in [0..NT-1] ] >;
// elliptic points
Ws := [ DiscMatrix(mrep, Quaternion(G ! q), CC ! p, CC) : q in G`ShimFDGens ];
pts := []; ords := [];
for W in Ws do
    tr := Real(W[1][1] + W[2][2]);
    if Abs(tr) ge 1.999 then continue; end if;
    ord := Abs(tr) lt 0.5 select 2 else 3;
    a := W[1][1]; b := W[1][2]; c := W[2][1]; d := W[2][2];
    dis := Sqrt((d-a)^2 + 4*b*c);
    for rt in [ (-(d-a)+dis)/(2*c), (-(d-a)-dis)/(2*c) ] do
        if Abs(rt) lt 0.999 then
            new := true; for q in pts do if Abs(q-rt) lt 1e-20 then new := false; end if; end for;
            if new then Append(~pts, rt); Append(~ords, ord); end if;
        end if;
    end for;
end for;
// Hauptmodul
P0 := CC ! (0.21 + 0.13*CC.1);
M1 := VerticalJoin(AutSystem(samp, 12, NT, CC), Matrix(CC, 1, NT, [ P0^n : n in [0..NT-1] ]));
scal := [ Sqrt(&+[ Abs(M1[i][j])^2 : i in [1..Nrows(M1)] ]) : j in [1..NT] ];
for j in [1..NT] do for i in [1..Nrows(M1)] do M1[i][j] := M1[i][j]/scal[j]; end for; end for;
MR := RealForm(M1); Sig, L, Rt := SingularValueDecomposition(MR); r := Min(Nrows(MR), Ncols(MR));
b1 := [ CC | (Rt[r][i]   + CC.1*Rt[r][NT+i])/scal[i]   : i in [1..NT] ];
b2 := [ CC | (Rt[r-2][i] + CC.1*Rt[r-2][NT+i])/scal[i] : i in [1..NT] ];
s := func< w | F(b1,w)/F(b2,w) >;
i3 := [i : i in [1..#pts] | ords[i] eq 3]; i2 := [i : i in [1..#pts] | ords[i] eq 2];
z1 := s(pts[i3[1]]); z2 := s(pts[i3[2]]); z3 := s(pts[i2[1]]);
// Mobius sending z1 -> 0, z2 -> oo, z3 -> i   (so that m matches Elkies' parametrisation)
mOf := func< w | CC.1 * ((s(w)-z1)/(s(w)-z2)) * ((z3-z2)/(z3-z1)) >;
for q in pts do
    if Abs(s(q)-z2) lt 1e-30 then printf "  m(elliptic) = infinity (as expected)\n";
    else printf "  m(elliptic) = %o\n", ComplexField(8)!mOf(q); end if;
end for;
// the disc -24 CM points: trace 0, reduced norm 6 in O  (Z[sqrt -6], disc -24).  Elkies: t = 0 there,
// i.e. m = +-1.  t = -((m^2-1)/(2m))^2.
B := Basis(O); found := [];
for c1,c2,c3,c4 in [-6..6] do
    mu := c1*B[1]+c2*B[2]+c3*B[3]+c4*B[4];
    if Trace(mu) eq 0 and Norm(mu) eq 6 then Append(~found, mu); end if;
    if #found ge 2 then break c1; end if;
end for;
if #found eq 0 then printf "no trace-0 norm-6 element found in O\n"; else
  // Elkies: the star Hauptmodul t vanishes at the CM points of discriminant -24 (m = +-1).
  // Magma's CMPoints hits an internal RingClassGroup error here, so take the fixed point of mu
  // directly: mu has trace 0 and reduced norm 6, so Z[mu] = Z[sqrt -6], discriminant -24.
  for mu in found do
    M := ChangeRing(mrep(mu), CC);
    a := M[1][1]; b := M[1][2]; c := M[2][1]; d := M[2][2];
    dis := Sqrt((d-a)^2 + 4*b*c);
    for zc in [ (-(d-a)+dis)/(2*c), (-(d-a)-dis)/(2*c) ] do
      if Imaginary(zc) le 0 then continue; end if;
      wc := (zc - CC!p)/(zc - ComplexConjugate(CC!p));
      red := ShimuraReduceUnit(G ! (O ! 1), Generators(G), DD : z0 := DD ! wc);
      del := red[1];
      Wd := DiscMatrix(mrep, Quaternion(del), CC ! p, CC);
      wr := (Wd[1][1]*wc + Wd[1][2])/(Wd[2][1]*wc + Wd[2][2]);
      mv := mOf(wr);
      tv := -((mv^2-1)/(2*mv))^2;
      printf "mu = %o : CM(disc -24) reduced to |w| = %o\n", mu, RealField(4)!Abs(wr);
      printf "   m = %o   (Elkies predicts m = +-1)\n", ComplexField(10)!mv;
      printf "   t = %o   (Elkies predicts t = 0)\n", ComplexField(10)!tv;
      printf "   |m|-1 = %o , |t| = %o\n", RealField(4)!(Abs(mv)-1), RealField(4)!Abs(tv);
    end for;
  end for;
end if;
printf "CMCHK DONE\n";
quit;
