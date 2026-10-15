load "/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/vw6/vw3.m";
SetColumns(0);
PREC := StringToInteger(PR); NT := StringToInteger(NTERMS); CC := ComplexField(PREC);
A := QuaternionAlgebra(6); O := QuaternionOrder(A, 1);
G := FuchsianGroup(O); DD := UnitDisc( : Precision := PREC);
FD := FundamentalDomain(G, DD);
rho := Max([Abs(ComplexValue(v)) : v in FD]);
H := UpperHalfPlane();
p := ComplexValue(DiscToPlane(H, DD ! 0) : Precision := PREC);
mrep := FuchsianMatrixRepresentation(A : Precision := PREC);
samp := ReducedSamples(G, O, DD, mrep, p, rho, 300, CC);
printf "rho %o, %o samples, N %o, prec %o\n", RealField(6)!rho, #samp, NT, PREC;

// ---- the four elliptic points: fixed points in the disc of the elliptic generators ----
Ws := [ DiscMatrix(mrep, Quaternion(G ! q), CC ! p, CC) : q in G`ShimFDGens ];
pts := []; ords := [];
for W in Ws do
    tr := Real(W[1][1] + W[2][2]);
    if Abs(tr) ge 1.999 then continue; end if;
    ord := Abs(tr) lt 0.5 select 2 else 3;
    a := W[1][1]; b := W[1][2]; c := W[2][1]; d := W[2][2];
    disc := Sqrt((d-a)^2 + 4*b*c);
    for rt in [ (-(d-a)+disc)/(2*c), (-(d-a)-disc)/(2*c) ] do
        if Abs(rt) lt 0.999 then
            new := true;
            for q in pts do if Abs(q - rt) lt 1e-20 then new := false; end if; end for;
            if new then Append(~pts, rt); Append(~ords, ord); end if;
        end if;
    end for;
end for;
F := func< bs, w | &+[ bs[n+1]*w^n : n in [0..NT-1] ] >;
printf "elliptic points found: %o\n", #pts;
for i in [1..#pts] do printf "  P%o order %o at w = %o |w| = %o\n", i, ords[i], ComplexField(8)!pts[i], RealField(5)!Abs(pts[i]); end for;

// ---- weight-12 forms, with an extra condition f(P0) = 0 so that a ratio has degree 1 ----
k := 12;
P0 := CC ! (0.21 + 0.13*CC.1);
M0 := AutSystem(samp, k, NT, CC);
extra := Matrix(CC, 1, NT, [ P0^n : n in [0..NT-1] ]);
M1 := VerticalJoin(M0, extra);
function NormaliseColumnsKeep(Mx)
    sc := [ Sqrt(&+[ Abs(Mx[i][j])^2 : i in [1..Nrows(Mx)] ]) : j in [1..Ncols(Mx)] ];
    for j in [1..Ncols(Mx)] do for i in [1..Nrows(Mx)] do Mx[i][j] := Mx[i][j]/sc[j]; end for; end for;
    return <Mx, sc>;
end function;
tmp := NormaliseColumnsKeep(M1); M2 := tmp[1]; scal := tmp[2];
sv := SingVals(M2);
printf "weight %o with one extra condition: sv[1..5] %o\n", k, [RealField(3)!s : s in sv[1..5]];
// two smallest right singular vectors
MR := RealForm(M2);
Sig, L, Rt := SingularValueDecomposition(MR);
r := Min(Nrows(MR), Ncols(MR));
function CoeffVec(row)
    return [ CC | (Rt[row][i] + CC.1*Rt[row][NT+i])/scal[i] : i in [1..NT] ];
end function;
b1 := CoeffVec(r); b2 := CoeffVec(r-2);     // singular values come in pairs for the complex matrix
printf "check both vanish at P0: %o %o\n", RealField(3)!Abs(F(b1,P0)), RealField(3)!Abs(F(b2,P0));
svals := [ F(b1, q)/F(b2, q) : q in pts ];
for i in [1..#pts] do printf "  s(P%o) [order %o] = %o\n", i, ords[i], ComplexField(10)!svals[i]; end for;
// cross-ratio with the two order-3 points first, then the two order-2 points
i3 := [i : i in [1..#pts] | ords[i] eq 3];
i2 := [i : i in [1..#pts] | ords[i] eq 2];
if #i3 eq 2 and #i2 eq 2 then
    z1 := svals[i3[1]]; z2 := svals[i3[2]]; z3 := svals[i2[1]]; z4 := svals[i2[2]];
    cr := ((z1-z3)*(z2-z4))/((z1-z4)*(z2-z3));
    printf "\nCROSS-RATIO (order3, order3 ; order2, order2) = %o\n", ComplexField(12)!cr;
    printf "Elkies prediction = -1 ; difference = %o\n", RealField(5)!Abs(cr - (-1));
else
    printf "unexpected elliptic-point classification: orders %o\n", ords;
end if;

// ---- valence check, NOT forced by symmetry ----
// The valence formula gives deg div(f_k) = (k/2)(2g-2+sum(1-1/e_i)) = k/6 here.  For k = 4 that is
// 2/3, which can only be carried by the two order-3 points (1/3 each), so the weight-4 form must
// vanish at both order-3 points and at neither order-2 point.
M4 := NormaliseColumnsKeep(AutSystem(samp, 4, NT, CC));
sv4 := SingVals(M4[1]);
printf "\nweight 4: sv[1..3] %o\n", [RealField(3)!x : x in sv4[1..3]];
MR4 := RealForm(M4[1]);
S4, L4, R4 := SingularValueDecomposition(MR4);
r4 := Min(Nrows(MR4), Ncols(MR4));
b4 := [ CC | (R4[r4][i] + CC.1*R4[r4][NT+i])/M4[2][i] : i in [1..NT] ];
nrm := Max([Abs(F(b4, q)) : q in pts]);
for i in [1..#pts] do
    printf "  |f4(P%o)| / max  = %o   [order %o]\n", i, RealField(4)!(Abs(F(b4,pts[i]))/nrm), ords[i];
end for;
printf "HAUPT DONE\n";
quit;
