// From dumped Nelson weight-2 values to the binary octic of a genus-3 hyperelliptic quotient.
//
//   magma -b IN:=ev35_c46.m hypfit.m < /dev/null
//
// The steps, and what each one would catch if it were wrong.
//
//  1. The conic through the canonical image.  Its kernel must be ONE-dimensional with a clean
//     gap; a ramp here means the points are clustered or the curve is not hyperelliptic.
//  2. A coordinate x on the conic.  Writing the conic's symmetric matrix as S^t S and then
//     changing to w2^2 = w1 w3 puts the canonical image in the form (1 : x : x^2), so x is the
//     degree-2 map to the projective line and w1 is the differential dx/y.
//  3. dx/dz by Cauchy's formula on each circle.  Computed twice, from all NP points and from
//     every other one; the two must agree, which is the statement that the circle is inside the
//     disc of analyticity and that the truncation of the trapezoidal rule is below the noise.
//     The mean of the circle values must also reproduce x at the centre, which is evaluated
//     independently and is not used in the derivative.
//  4. y = (dx/dz)/w1 is anti-invariant under the hyperelliptic involution, so y^2 is a function
//     of x, and of degree 8.  Nine of the centres determine the octic and EVERY FURTHER CENTRE
//     IS A PREDICTION -- this is the step that would fail outright if the curve were not
//     hyperelliptic, or if the conic were an artefact of clustering.
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/nelson/wk.m";
SetColumns(0);
PRLOAD := StringToInteger(assigned PW select PW else "0");

eval Read(IN);
prec := (PRLOAD gt 0) select PRLOAD else PRECV;
CC := ComplexField(prec); RR := RealField(prec); Ic := CC.1;
R4 := RealField(6);
g := GG;
error if g ne 3, "this script does genus 3";

np := NPC; nc := NCIRC; rho := RR ! RHO;
pts := [ CC | PTS[j][1] + Ic*PTS[j][2] : j in [1..#PTS] ];
// the holomorphic value: Theorem 3.1 returns Im(z1) Im(z) conj(f(z1)) f(z), so dividing by Im(z)
// leaves a constant multiple of f(z), and the constant is the same at every z
Gv := [ [ CC | (VALS[j][i][1] + Ic*VALS[j][i][2]) / Imaginary(pts[j]) : i in [1..g] ]
        : j in [1..#PTS] ];
printf "=== X_0^%o(%o)/%o  genus %o, cut %o, %o circles of %o, rho %o, working at %o digits\n",
       DD, NN, WW, g, CUTV, nc, np, R4!rho, prec;
printf "Atkin-Lehner eigenvalues of the %o differentials at %o:\n", g, HALLS;
for r in ALSIG do printf "   %o\n", r; end for;
hypAL := [ HALLS[i] : i in [1..#HALLS] | forall{ r : r in ALSIG | r[i] eq -1 } and HALLS[i] ne 1 ];
printf "Hall divisors acting as -1 on every differential (an Atkin-Lehner hyperelliptic "
     * "involution would be one of these): %o\n", hypAL;

// ---------- 1. the conic ----------
// use a well-spread subset: the centres and one point of each circle
idx := [ j : j in [1..#PTS] | ISCTR[j] eq 1 ] cat
       [ j : j in [1..#PTS] | ISCTR[j] eq 0 and (j mod np) eq 3 ];
mons := [ [a,b] : b in [a..g], a in [1..g] ];
M := Matrix(CC, #idx, #mons,
     &cat[ [ Gv[j][m[1]]*Gv[j][m[2]] : m in mons ] : j in idx ]);
for i in [1..#idx] do
    nm := Sqrt(&+[ Abs(M[i][c])^2 : c in [1..#mons] ]);
    for c in [1..#mons] do M[i][c] := M[i][c]/nm; end for;
end for;
for c in [1..#mons] do
    nm := Sqrt(&+[ Abs(M[i][c])^2 : i in [1..#idx] ]);
    for i in [1..#idx] do M[i][c] := M[i][c]/nm; end for;
end for;
n2 := #idx;
MR := Matrix(RR, 2*n2, 2*#mons, [ RR ! 0 : i in [1..4*n2*#mons] ]);
for i in [1..n2] do for j in [1..#mons] do
    a := Real(M[i][j]); b := Imaginary(M[i][j]);
    MR[i][j] := a; MR[i][#mons+j] := -b; MR[n2+i][j] := b; MR[n2+i][#mons+j] := a;
end for; end for;
Sg, U, V := SingularValueDecomposition(MR);
r := Min(Nrows(MR), Ncols(MR));
sv := Reverse([ Abs(Sg[i][i]) : i in [1..r] ]);
sv := [ sv[i] : i in [1..#sv by 2] ];
printf "conic: singular values %o\n", [ R4!x : x in sv ];
printf "       residual %o, next %o, ratio %o\n", R4!sv[1], R4!sv[2], R4!(sv[2]/sv[1]);
error if sv[2]/sv[1] lt 10^4,
   "no clean gap: either the points are clustered or there is no conic (not hyperelliptic)";
// the kernel vector is the last column of V (real form: real part then imaginary part)
kcol := Ncols(V);
kv := [ CC | V[j][kcol] + Ic*V[#mons+j][kcol] : j in [1..#mons] ];
// undo the column normalisation applied above
cnm := [ Sqrt(&+[ Abs(Matrix(CC,1,1,[Gv[j][mons[c][1]]*Gv[j][mons[c][2]]])[1][1])^2 : j in idx ])
         : c in [1..#mons] ];
// (recompute the column norms on the row-normalised matrix exactly as used)
Mx := Matrix(CC, #idx, #mons,
     &cat[ [ Gv[j][m[1]]*Gv[j][m[2]] : m in mons ] : j in idx ]);
for i in [1..#idx] do
    nm := Sqrt(&+[ Abs(Mx[i][c])^2 : c in [1..#mons] ]);
    for c in [1..#mons] do Mx[i][c] := Mx[i][c]/nm; end for;
end for;
cn := [ Sqrt(&+[ Abs(Mx[i][c])^2 : i in [1..#idx] ]) : c in [1..#mons] ];
qc := [ CC | kv[c]/cn[c] : c in [1..#mons] ];
// symmetric matrix of the quadratic form
Q := ZeroMatrix(CC, 3, 3);
for c in [1..#mons] do
    a := mons[c][1]; b := mons[c][2];
    if a eq b then Q[a][a] := qc[c];
    else Q[a][b] := qc[c]/2; Q[b][a] := qc[c]/2; end if;
end for;
// how well the conic holds at the centres, as a relative figure
rels := [ RR | ];
for j in [1..#PTS] do
    v := Matrix(CC,3,1,[Gv[j][1],Gv[j][2],Gv[j][3]]);
    num := Abs((Transpose(v)*Q*v)[1][1]);
    den := &+[ Abs(Q[a][b])*Abs(Gv[j][a])*Abs(Gv[j][b]) : a,b in [1..3] ];
    Append(~rels, num/den);
end for;
printf "conic relative residual over all %o points: max %o, median %o\n",
       #PTS, R4!Maximum(rels), R4!Sort(rels)[#rels div 2];

// ---------- 2. a coordinate on the conic ----------
// Q = S^t S by symmetric (LDL^t) reduction, then w1 = i u1 - u3, w2 = u2, w3 = i u1 + u3
// turns u1^2+u2^2+u3^2 into w2^2 - w1 w3, so the canonical image is (w1:w2:w3) = (1:x:x^2).
A := Q; L := IdentityMatrix(CC, 3); Dg := [ CC | 0,0,0 ];
B := A;
for i in [1..3] do
    Dg[i] := B[i][i];
    error if Abs(Dg[i]) lt 10^(-prec div 2), "the conic is degenerate at the symmetric reduction";
    for j in [i+1..3] do
        f := B[j][i]/Dg[i];
        L[j][i] := f;
        for c in [1..3] do B[j][c] := B[j][c] - f*B[i][c]; end for;
    end for;
    for j in [i+1..3] do
        f := L[j][i];
        for c in [1..3] do B[c][j] := B[c][j] - f*B[c][i]; end for;
    end for;
end for;
S := DiagonalMatrix(CC, [ Sqrt(Dg[i]) : i in [1..3] ]) * Transpose(L);
chk := Transpose(S)*S - Q;
printf "symmetric factorisation residual: %o\n",
       R4!(Maximum([ Abs(chk[a][b]) : a,b in [1..3] ])/Maximum([ Abs(Q[a][b]) : a,b in [1..3] ]));
Amat := Matrix(CC, 3, 3, [ Ic,0,-1, 0,1,0, Ic,0,1 ]);
T := Amat*S;
function Wof(j)
    v := Matrix(CC,3,1,[Gv[j][1],Gv[j][2],Gv[j][3]]);
    u := T*v;
    return u[1][1], u[2][1], u[3][1];
end function;
// the normal form must hold: w2^2 = w1 w3
nf := [ RR | ];
for j in [1..#PTS] do
    w1,w2,w3 := Wof(j);
    Append(~nf, Abs(w2^2 - w1*w3)/(Abs(w2)^2 + Abs(w1*w3)));
end for;
printf "normal form w2^2 = w1 w3, relative residual: max %o\n", R4!Maximum(nf);

// ---------- 3. x and dx/dz ----------
xs := [ CC | ]; ys := [ CC | ]; dchk := [ RR | ]; cchk := [ RR | ];
for c in [1..nc] do
    base := (c-1)*(np+1) + 1;            // the centre
    cir := [ base + t : t in [1..np] ];  // the circle, t = 0 .. np-1 in order
    xc := [ CC | ];
    for j in cir do
        w1,w2,w3 := Wof(j);
        Append(~xc, w2/w1);
    end for;
    w1c, w2c, w3c := Wof(base);
    xcen := w2c/w1c;
    om := Exp(2*Pi(CC)*Ic/np);
    a1 := (&+[ xc[t+1]*om^(-t) : t in [0..np-1] ])/(np*rho);
    a0 := (&+[ xc[t+1] : t in [0..np-1] ])/np;
    // the same derivative from half the points: a genuine convergence check on the circle
    a1h := (&+[ xc[2*t+1]*om^(-2*t) : t in [0..(np div 2)-1] ])/((np div 2)*rho);
    Append(~dchk, Abs(a1-a1h)/Abs(a1));
    Append(~cchk, Abs(a0-xcen)/Abs(xcen));
    Append(~xs, xcen);
    Append(~ys, a1/w1c);
end for;
printf "Cauchy derivative, %o points against %o: max relative difference %o\n",
       np, np div 2, R4!Maximum(dchk);
printf "circle mean against the independently evaluated centre: max relative difference %o\n",
       R4!Maximum(cchk);

// ---------- 4. the octic ----------
// recentre and rescale x, which is an affine change of the coordinate and sends y to y/sigma,
// so the shape y^2 = (degree 8) is kept and the Vandermonde system is conditioned
mu := (&+xs)/nc;
sig := Sqrt((&+[ Abs(x-mu)^2 : x in xs ])/nc);
xt := [ (x-mu)/sig : x in xs ];
yt := [ y/sig : y in ys ];
printf "x after recentring: moduli %o\n", [ R4!Abs(x) : x in xt ];

nun := 9;
error if nc lt nun+1, "not enough centres to leave a prediction";
Amx := Matrix(CC, nun, nun, &cat[ [ xt[j]^(i-1) : i in [1..nun] ] : j in [1..nun] ]);
bv := Matrix(CC, nun, 1, [ yt[j]^2 : j in [1..nun] ]);
cf := Amx^(-1)*bv;
coef := [ CC | cf[i][1] : i in [1..nun] ];
printf "octic fitted from %o centres; %o further centres are PREDICTIONS:\n", nun, nc-nun;
maxpr := RR ! 0;
for j in [nun+1..nc] do
    pv := &+[ coef[i]*xt[j]^(i-1) : i in [1..nun] ];
    e := Abs(pv - yt[j]^2)/Abs(yt[j]^2);
    maxpr := Maximum(maxpr, e);
    printf "   centre %o: predicted/observed relative difference %o\n", j, R4!e;
end for;
printf "WORST PREDICTION: %o\n", R4!maxpr;
printf "octic coefficients (should be real up to the noise if the curve is defined over R "
     * "in this coordinate -- they need not be):\n";
for i in [1..nun] do printf "   c%o = %o\n", i-1, coef[i]; end for;

PC<Xv> := PolynomialRing(CC);
foct := &+[ coef[i]*Xv^(i-1) : i in [1..nun] ];
printf "leading coefficient %o, so the degree is %o\n", R4!Abs(coef[9]),
       Abs(coef[9]) lt 10^(-prec div 3) select 7 else 8;
rts := Roots(foct);
printf "roots of the octic (the branch points), %o of them with multiplicities %o\n",
       #rts, [ t[2] : t in rts ];
for t in rts do printf "   %o\n", t[1]; end for;

// Shioda invariants are the coordinate-free data of a binary octic, so they are what can be
// recognised as rationals; they are printed here and recognised in the next script.
try
    SI := ShiodaInvariants(foct);
    printf "Shioda invariants:\n";
    for i in [1..#SI] do printf "   J%o = %o\n", i+1, SI[i]; end for;
catch e
    printf "ShiodaInvariants failed: %o\n", e`Object;
end try;
printf "HYPFIT DONE\n";
quit;
