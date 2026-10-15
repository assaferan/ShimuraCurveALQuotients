// Is the disc automorphy factor used before equal to the paper's half-plane one with the (1-w)^k
// factors?  Paper eq (9): (1-w')^k sum b_n w'^n = j(g,z)^k (1-w)^k sum b_n w^n, j the HALF-PLANE
// factor cz+d.  Earlier code used F(w') = j_D(W,w)^k F(w) with F = sum a_n w^n.
load "/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/vw/vwlib.m";
SetColumns(0);
PREC := 60; CC := ComplexField(PREC);
A := QuaternionAlgebra(6); O := QuaternionOrder(A, 1);
G := FuchsianGroup(O);
DD := UnitDisc( : Precision := PREC);
FD := FundamentalDomain(G, DD);
rho := Max([Abs(ComplexValue(v)) : v in FD]);
H := UpperHalfPlane();
p := ComplexValue(DiscToPlane(H, DD ! 0) : Precision := PREC);
mrep := FuchsianMatrixRepresentation(A : Precision := PREC);
printf "rho = %o, p = %o\n", RealField(6)!rho, ComplexField(6)!p;
zOf := func< w | (ComplexConjugate(CC!p)*w - CC!p)/(w - 1) >;
wOf := func< z | (z - CC!p)/(z - ComplexConjugate(CC!p)) >;
for gi in [1, 2] do
    q := Quaternion(G ! G`ShimFDGens[gi]);
    M := ChangeRing(mrep(q), CC);
    W := DiscMatrix(mrep, q, CC!p, CC);
    printf "\n--- generator %o: det(M)=%o det(W)=%o\n", gi, ComplexField(6)!Determinant(M), ComplexField(6)!Determinant(W);
    for t in [0.3, 0.45] do
        w := CC ! (t * Exp(2*Pi(RealField(PREC))*CC.1*0.37));
        z := zOf(w);
        gz := (M[1][1]*z + M[1][2])/(M[2][1]*z + M[2][2]);
        wp := DiscAct(W, w);
        printf "  |w|=%o  w(gz) vs DiscAct agree to %o\n", RealField(4)!Abs(w), RealField(3)!Abs(wOf(gz) - wp);
        jH := M[2][1]*z + M[2][2];
        jD := DiscJ(W, w);
        for k in [4, 6] do
            lhs := jD^k;                                   // earlier code's factor
            rhs := jH^k * (1-w)^k / (1-wp)^k;              // paper's effective factor
            printf "    k=%o: jD^k/[jH^k(1-w)^k/(1-w')^k] = %o\n", k, ComplexField(8)!(lhs/rhs);
        end for;
    end for;
end for;
printf "DBGB DONE\n";
quit;
