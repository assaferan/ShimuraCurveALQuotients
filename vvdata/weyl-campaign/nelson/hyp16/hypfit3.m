// Nelson values -> the equation of a genus-3 hyperelliptic Atkin-Lehner quotient.
//
//   magma -b IN:=ev35_c46.m hypfit2.m < /dev/null
//
// THE COORDINATE IS PINNED BY ATKIN-LEHNER, WHICH IS WHAT MAKES THIS RECOGNISABLE.
// The hyperelliptic involution is central in the automorphism group, so an Atkin-Lehner
// involution u descends to an involution of the projective line X/iota and the coordinate can be
// chosen with u : x -> -x.  The octic is then EVEN, f(x) = g(x^2), and X/u is the genus-1 curve
// Y^2 = g(t).  Both of these are checks and not assumptions: the evenness is read off the fit,
// and the quadric must have no cross terms between the u-invariant differential and the other
// two, which is forced by u-equivariance and would fail if the eigenvalue convention were wrong.
//
// The point of the normalisation is HEIGHT.  The coordinate-free Shioda invariants of a binary
// octic carry sixty-digit numerators at the control, which no numerical evaluation here could
// recognise.  The quartic g with its two marked points t = 0, infinity has instead the three
// invariants c0c2/c1^2, c1c3/c2^2, c2c4/c3^2, unchanged by t -> lambda t and by rescaling g; at
// X_0^35(1) they are 1939/952200, 34500/76729 and 1939/5000, which fifteen digits settle.
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/nelson/wk.m";
SetColumns(0);
dmp := eval (Read(IN) cat
   "\nreturn <DD, NN, WW, GG, PRECV, CUTV, RHO, NCIRC, NPC, HALLS, ALSIG, PTS, OWNER, ISCTR, VALS>;");
DD := dmp[1]; NN := dmp[2]; WW := dmp[3]; GG := dmp[4]; PRECV := dmp[5]; CUTV := dmp[6];
RHO := dmp[7]; NCIRC := dmp[8]; NPC := dmp[9]; HALLS := dmp[10]; ALSIG := dmp[11];
PTS := dmp[12]; OWNER := dmp[13]; ISCTR := dmp[14]; VALS := dmp[15];
prec := PRECV;
CC := ComplexField(prec); RR := RealField(prec); Ic := CC.1;
R6 := RealField(8);
g := GG;
error if g ne 3, "genus 3 only in this script";
np := NPC; nc := NCIRC; rho := RR ! RHO;
pts := [ CC | PTS[j][1] + Ic*PTS[j][2] : j in [1..#PTS] ];
Gv := [ [ CC | (VALS[j][i][1] + Ic*VALS[j][i][2]) / Imaginary(pts[j]) : i in [1..g] ]
        : j in [1..#PTS] ];
printf "=== X_0^%o(%o)/%o  genus %o, cut %o, %o digits, %o circles of %o, rho %o\n",
       DD, NN, WW, g, CUTV, prec, nc, np, R6!rho;

// ---------- the Atkin-Lehner involution with a genus-1 quotient ----------
// hypfit3.m = hypfit2.m plus the optional argument UM, which names the Atkin-Lehner involution
// u to pin the coordinate by when several have a genus-1 quotient (X_0^1326(1) has three).
// Without UM the choice is hypfit2.m's: the first Hall divisor with exactly one +1 form.
uix := 0;
if assigned UM then
    uix := Index(HALLS, StringToInteger(UM));
    error if uix eq 0, "UM is not a Hall divisor";
    error if #[ r : r in ALSIG | r[uix] eq 1 ] ne 1, "UM does not have a genus-1 quotient (needs exactly one +1 form)";
else
for i in [1..#HALLS] do
    if HALLS[i] eq 1 then continue; end if;
    np1 := #[ r : r in ALSIG | r[i] eq 1 ];
    if np1 eq 1 then uix := i; break; end if;
end for;
end if;
error if uix eq 0, "no Atkin-Lehner involution with a genus-1 quotient; this script needs one";
um := HALLS[uix];
pind := [ j : j in [1..g] | ALSIG[j][uix] eq 1 ][1];
aind := [ j : j in [1..g] | ALSIG[j][uix] eq -1 ];
printf "Atkin-Lehner involution w_%o : quotient genus 1; invariant differential #%o, "
     * "anti-invariant #%o\n", um, pind, aind;

// ---------- the conic ----------
idx := [ j : j in [1..#PTS] | ISCTR[j] eq 1 ] cat
       [ j : j in [1..#PTS] | ISCTR[j] eq 0 and (j mod np) eq 3 ];
mons := [ [a,b] : b in [a..g], a in [1..g] ];
Mx := Matrix(CC, #idx, #mons, &cat[ [ Gv[j][m[1]]*Gv[j][m[2]] : m in mons ] : j in idx ]);
for i in [1..#idx] do
    nm := Sqrt(&+[ Abs(Mx[i][c])^2 : c in [1..#mons] ]);
    for c in [1..#mons] do Mx[i][c] := Mx[i][c]/nm; end for;
end for;
cn := [ Sqrt(&+[ Abs(Mx[i][c])^2 : i in [1..#idx] ]) : c in [1..#mons] ];
My := Mx;
for c in [1..#mons] do for i in [1..#idx] do My[i][c] := My[i][c]/cn[c]; end for; end for;
n2 := #idx;
MR := Matrix(RR, 2*n2, 2*#mons, [ RR ! 0 : i in [1..4*n2*#mons] ]);
for i in [1..n2] do for j in [1..#mons] do
    a := Real(My[i][j]); b := Imaginary(My[i][j]);
    MR[i][j] := a; MR[i][#mons+j] := -b; MR[n2+i][j] := b; MR[n2+i][#mons+j] := a;
end for; end for;
Sg, U, V := SingularValueDecomposition(MR);
rr := Min(Nrows(MR), Ncols(MR));
sv := Reverse([ Abs(Sg[i][i]) : i in [1..rr] ]);
sv := [ sv[i] : i in [1..#sv by 2] ];
printf "conic singular values: %o\n", [ R6!x : x in sv ];
printf "conic residual %o, next %o, ratio %o\n", R6!sv[1], R6!sv[2], R6!(sv[2]/sv[1]);
error if sv[2]/sv[1] lt 10^4, "no clean conic";
// Magma's documented convention for which of the two orthogonal factors holds the right
// singular vectors, and whether they are its rows or its columns, is not something to guess:
// both readings are tried and the one that actually annihilates the matrix is taken.
kcol := Ncols(V);
function KernRes(vec)
    s := RR ! 0;
    for i in [1..#idx] do
        s +:= Abs(&+[ My[i][c]*vec[c] : c in [1..#mons] ])^2;
    end for;
    return Sqrt(s)/Sqrt(&+[ Abs(c)^2 : c in vec ]);
end function;
vcol := [ CC | V[j][kcol] + Ic*V[#mons+j][kcol] : j in [1..#mons] ];
vrow := [ CC | V[kcol][j] + Ic*V[kcol][#mons+j] : j in [1..#mons] ];
rc := KernRes(vcol); rw := KernRes(vrow);
printf "kernel vector: residual as a column of V %o, as a row %o -- taking the %o\n",
       R6!rc, R6!rw, rc lt rw select "column" else "row";
kvec := rc lt rw select vcol else vrow;
error if Minimum(rc,rw) gt 10^(-8), "neither reading of the SVD gives a kernel vector";
qc := [ CC | kvec[j]/cn[j] : j in [1..#mons] ];
Q := ZeroMatrix(CC, 3, 3);
for c in [1..#mons] do
    a := mons[c][1]; b := mons[c][2];
    if a eq b then Q[a][a] := qc[c]; else Q[a][b] := qc[c]/2; Q[b][a] := qc[c]/2; end if;
end for;
scl := Maximum([ Abs(Q[a][b]) : a,b in [1..3] ]);
// u-equivariance: the quadric cannot pair the invariant differential with an anti-invariant one
printf "u-equivariance of the quadric (these must vanish): %o  %o\n",
       R6!(Abs(Q[pind][aind[1]])/scl), R6!(Abs(Q[pind][aind[2]])/scl);

// ---------- the coordinate ----------
// Q = qpp ep^2 + (binary form in the two anti-invariant coordinates); the conic says
// ep^2 = (-qaa/qpp)(ea - r1 eb)(ea - r2 eb), so w2 = ep, w1 = ea - r1 eb, w3 = that multiple.
a1 := aind[1]; a2 := aind[2];
qpp := Q[pind][pind]; qaa := Q[a1][a1]; qab := Q[a1][a2]; qbb := Q[a2][a2];
disc := qab^2 - qaa*qbb;
r1 := (-qab + Sqrt(disc))/qaa; r2 := (-qab - Sqrt(disc))/qaa;
kap := -qaa/qpp;
function Wof(j)
    ep := Gv[j][pind]; ea := Gv[j][a1]; eb := Gv[j][a2];
    w1 := ea - r1*eb;
    w3 := kap*(ea - r2*eb);
    return w1, ep, w3;
end function;
nf := [ RR | ];
for j in [1..#PTS] do
    w1,w2,w3 := Wof(j);
    Append(~nf, Abs(w2^2 - w1*w3)/(Abs(w2)^2 + Abs(w1*w3)));
end for;
printf "normal form w2^2 = w1 w3: max relative residual %o\n", R6!Maximum(nf);

// ---------- x, dx/dz, y ----------
xs := [ CC | ]; ys := [ CC | ]; dchk := [ RR | ]; cchk := [ RR | ];
for c in [1..nc] do
    base := (c-1)*(np+1) + 1;
    xc := [ CC | ];
    for t in [1..np] do
        w1,w2,w3 := Wof(base+t);
        Append(~xc, w2/w1);
    end for;
    w1c,w2c,w3c := Wof(base);
    xcen := w2c/w1c;
    om := Exp(2*Pi(CC)*Ic/np);
    d1 := (&+[ xc[t+1]*om^(-t) : t in [0..np-1] ])/(np*rho);
    d0 := (&+[ xc[t+1] : t in [0..np-1] ])/np;
    d1h := (&+[ xc[2*t+1]*om^(-2*t) : t in [0..(np div 2)-1] ])/((np div 2)*rho);
    Append(~dchk, Abs(d1-d1h)/Abs(d1));
    Append(~cchk, Abs(d0-xcen)/Abs(xcen));
    Append(~xs, xcen); Append(~ys, d1/w1c);
end for;
printf "Cauchy derivative %o points against %o: max relative difference %o\n",
       np, np div 2, R6!Maximum(dchk);
printf "circle mean against independently evaluated centre: max relative difference %o\n",
       R6!Maximum(cchk);
printf "x at the centres: %o\n", [ R6!Abs(x) : x in xs ];
// The overall scale of x is the one freedom Atkin-Lehner leaves (x -> lambda x keeps u : x -> -x),
// so rescaling to modulus about 1 is free and it is what keeps the degree-8 Vandermonde system
// conditioned.  The three invariants below are unchanged by it: c_i picks up sigma^(2i) and they
// are the combinations in which that cancels.
sigx := Exp((&+[ Log(Abs(xx)) : xx in xs ])/nc);
printf "rescaling x by %o for conditioning (A, B, C are unchanged by this)\n", R6!sigx;
xs := [ xx/sigx : xx in xs ];
ys := [ yy/sigx : yy in ys ];

// ---------- the octic ----------
// no recentring: the coordinate is pinned by Atkin-Lehner (u : x -> -x), only x -> lambda x
// and x -> 1/x remain, and the octic must come out EVEN.
// The octic is even, so y^2 is a QUARTIC in t = x^2: five unknowns rather than nine, which is a
// far better conditioned Vandermonde and leaves nine predictions instead of five.  That the
// quartic fit succeeds at all IS the statement that the octic is even, hence that u acts as
// x -> -x; a curve whose octic had odd terms could not satisfy it at nine further centres.
ts := [ xx^2 : xx in xs ];
nun := 5;
error if nc lt nun+1, "not enough centres to leave a prediction";
Amx := Matrix(CC, nun, nun, &cat[ [ ts[j]^(i-1) : i in [1..nun] ] : j in [1..nun] ]);
bv := Matrix(CC, nun, 1, [ ys[j]^2 : j in [1..nun] ]);
cf := Amx^(-1)*bv;
coef := [ CC | cf[i][1] : i in [1..nun] ];
maxpr := RR ! 0;
printf "quartic in t = x^2 fitted on %o centres, %o PREDICTIONS:\n", nun, nc-nun;
for j in [nun+1..nc] do
    pv := &+[ coef[i]*ts[j]^(i-1) : i in [1..nun] ];
    e := Abs(pv - ys[j]^2)/Abs(ys[j]^2);
    maxpr := Maximum(maxpr, e);
    printf "   centre %o: relative difference %o\n", j, R6!e;
end for;
printf "WORST PREDICTION (quartic in x^2): %o\n", R6!maxpr;

// the full octic as a diagnostic: its odd coefficients measure the same thing directly
if nc ge 10 then
    Ao := Matrix(CC, 9, 9, &cat[ [ xs[j]^(i-1) : i in [1..9] ] : j in [1..9] ]);
    bo := Matrix(CC, 9, 1, [ ys[j]^2 : j in [1..9] ]);
    co := Ao^(-1)*bo;
    oc := [ CC | co[i][1] : i in [1..9] ];
    sco := Maximum([ Abs(c) : c in oc ]);
    printf "octic fitted independently on 9 centres; odd coefficients relative to the largest "
         * "(these measure the evenness directly): %o\n",
           [ R6!(Abs(oc[i])/sco) : i in [2,4,6,8] ];
end if;

// The five-point solve above is the honest prediction test, but it throws away nine centres and
// the Vandermonde in t is ill-conditioned because |t| spans two or three orders of magnitude.
// For the VALUES, use every centre, each equation divided by y^2 so that all residuals are
// relative and no centre dominates merely by being large.
Aw := Matrix(CC, nc, nun, &cat[ [ ts[j]^(i-1)/ys[j]^2 : i in [1..nun] ] : j in [1..nc] ]);
bw := Matrix(CC, nc, 1, [ CC | 1 : j in [1..nc] ]);
Ah := Matrix(CC, nun, nc, &cat[ [ ComplexConjugate(Aw[j][i]) : j in [1..nc] ] : i in [1..nun] ]);
cls := (Ah*Aw)^(-1)*(Ah*bw);
coefL := [ CC | cls[i][1] : i in [1..nun] ];
resL := [ RR | Abs(&+[ coefL[i]*ts[j]^(i-1) : i in [1..nun] ]/ys[j]^2 - 1) : j in [1..nc] ];
printf "weighted least squares over all %o centres: max relative residual %o\n",
       nc, R6!Maximum(resL);
printf "five-point and least-squares quartics agree to: %o\n",
       R6!Maximum([ Abs(coefL[i]-coef[i])/Abs(coef[i]) : i in [1..nun] ]);

// ---------- the marked invariants of the quartic ----------
gc := coefL;
A := gc[1]*gc[3]/gc[2]^2; B := gc[2]*gc[4]/gc[3]^2; Cc := gc[3]*gc[5]/gc[4]^2;
printf "\nmarked-quartic invariants (these are the rational numbers):\n";
printf "   A = %o\n   B = %o\n   C = %o\n", A, B, Cc;
printf "   imaginary parts (a measure of the noise): %o  %o  %o\n",
       R6!Abs(Imaginary(A)), R6!Abs(Imaginary(B)), R6!Abs(Imaginary(Cc));

function Recognise(z, acc)
    r := Real(z);
    scale := Abs(r) gt 1 select Abs(r) else Parent(r) ! 1;
    eps := Max(10^(-acc) * scale, 20 * Abs(Imaginary(z)));
    for e in [2..acc] do
        q := BestApproximation(r, 10^e);
        if Abs(r - q) lt eps and Denominator(q)^2 * eps lt 1/2 then return true, q; end if;
    end for;
    return false, 0;
end function;
ACC := StringToInteger(assigned ACCD select ACCD else Sprint(Maximum(8, prec div 3)));
printf "\nrecognising with the classical certainty criterion at %o digits:\n", ACC;
for nm in [ <"A",A>, <"B",B>, <"C",Cc> ] do
    ok, q := Recognise(nm[2], ACC);
    printf "   %o = %o\n", nm[1], ok select Sprint(q) else "NOT RECOGNISED";
end for;
// swapping the two roots of the binary form is x -> 1/x, which reverses the quartic and so
// exchanges A and C; both orderings are legitimate and the arithmetic stage can take either
printf "the other root ordering (x -> 1/x) gives (A,B,C) = (C,B,A)\n";
printf "HYPFIT3 DONE\n";
quit;
