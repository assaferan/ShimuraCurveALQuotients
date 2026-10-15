// Diagnose the four exact null vectors at k = 0 mod 4.  Build the prior system, find its
// nullspace, then test each null vector's automorphy at FRESH points not in the sample.
load "/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/vw/vwlib.m";
SetColumns(0);
PREC := 80; CC := ComplexField(PREC); NT := 40;
A := QuaternionAlgebra(6); O := QuaternionOrder(A, 1);
G := FuchsianGroup(O); sig := Signature(G);
printf "signature %o\n", sig;
DD := UnitDisc( : Precision := PREC);
FD := FundamentalDomain(G, DD);
R := Max([Abs(ComplexValue(v)) : v in FD]) * 1.02;
H := UpperHalfPlane();
p := ComplexValue(DiscToPlane(H, DD ! 0) : Precision := PREC);
mrep := FuchsianMatrixRepresentation(A : Precision := PREC);
gens := [Quaternion(G ! q) : q in G`ShimFDGens];
Ws := [ DiscMatrix(mrep, q, CC ! p, CC) : q in gens ];
printf "R = %.5o, %o generators\n", R, #Ws;
// orders of the generators in the group, and the SU(1,1) check
for i->W in Ws do
    tr := W[1][1] + W[2][2];
    printf "  gen %o: trace %o, |c| %o\n", i, ComplexField(6)!tr, RealField(5)!Abs(W[2][1]);
end for;
for k in [4, 6] do
    Mx0, maxm, used := AutomorphySystem(Ws, k, NT, 400, R, CC);
    printf "\nk=%o: %o rows, %o cols, maxmod %.4o\n", k, Nrows(Mx0), NT, maxm;
    // how many rows are (numerically) zero?
    nz := 0;
    for i in [1..Nrows(Mx0)] do
        if &+[Abs(Mx0[i][j]) : j in [1..NT]] lt 1e-60 then nz +:= 1; end if;
    end for;
    printf "  identically-zero rows: %o of %o\n", nz, Nrows(Mx0);
    Mx := NormaliseColumns(Mx0);
    sv := SingVals(Mx);
    printf "  smallest 8 singular values: %o\n", [RealField(3)!s : s in sv[1..Min(8,#sv)]];
end for;
printf "DBGA DONE\n";
quit;
