// Stage 3 at squarefree level N > 1 with the OLDFORM part of the space supplied.
//
// Theorem 3.1 evaluates a form on X_0^D(M) from a newform of level exactly D*M.  The
// differentials of X_0^D(N) are the D-new part of S_2(Gamma_0(DN)), which is the direct sum over
// M | N of the oldclasses of the newforms g of level D*M; the oldclass of g is spanned by the
// pullbacks g|w_S, S running over the subsets of the primes T dividing N/M, where w_S is an
// element of the level-N Eichler order O_N of reduced norm prod(S) normalising O_N (a product of
// the Atkin-Lehner elements w_p, p in S) and f|w is the weight-2 slash  (f|w)(z) = nrd(w)
// j(w,z)^-2 f(wz), which is the pullback of the differential f dz.  On that span w_q for q | DM
// acts by the quaternionic eigenvalue eps_q(g) (a_g(q) at q | D, -a_g(q) at q | M) and w_p for p
// in T permutes S -> S xor {p} with no scalar, because w_p^2 is p times a unit of norm one; so the
// eigenvectors are v_chi = sum_S chi(S) g|w_S for the characters chi of {+-1}^T, with w_m
// eigenvalue prod_{q | m, q | DM} eps_q(g) prod_{p | m, p in T} chi(p).
//
// What the evaluator supplies at level M is S_g(z1, z2) = y1 y2 conj(g_B(z1)) g_B(z2), and the
// same product for the pullback is, with j = j(w_S, z2),
//
//     y1 y2 conj(g_B(z1)) (g_B|w_S)(z2) = (conj(j)/j) S_g(z1, w_S z2),
//
// so every member of the oldclass is read off the level-M lattice at the moved point, with the
// same constant conj(g_B(z1)) in front as the untransformed pullback.  No further normalisation
// enters: the relative scale of the pullbacks is fixed by the slash operator, and that is what
// makes v_chi an eigenvector with eigenvalue exactly chi(p).
//
// The orders: O_N = QuaternionOrder(A, N) as in nelson3.m, O_1 a maximal order containing it,
// O_M = O_N + M O_1, which is the Eichler order of level M between them; w_m is found in the
// two-sided ideal J_m = {x in O_N : trd(x y) in m Z for all y in O_N} (the radical of O_N/mO_N)
// as an element of reduced norm exactly m, which forces it to generate J_m and so to normalise.
//
// Checks: the W-fixed count including oldforms must equal the genus (eigenvalue bookkeeping), and
// with CHECKW:=n the W-invariance of every column is measured at n points, V(w_m z) conj(j)/j
// against V(z) for each generator m of W, with the rejected characters of each oldclass carried
// along so that their predicted eigenvalues (not 1) are measured too.
//
//   magma -b DISC:=46 LEV:=3 WS:=1,3 PR:=30 CUT:=40 NPTS:=18 CHECKW:=2 stage3o.m < /dev/null
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/nelson/nelson3.m";

D := StringToInteger(DISC);
N := (assigned LEV) select StringToInteger(LEV) else 1;
W := { StringToInteger(x) : x in Split(WS, ",") };
PREC := StringToInteger(PR); cutv := StringToInteger(CUT); npts := StringToInteger(NPTS);
nchk := (assigned CHECKW) select StringToInteger(CHECKW) else 0;
CC := ComplexField(PREC); RR := RealField(PREC); Ic := CC.1;
C8 := ComplexField(8);
k := 2;

g := GenusShimuraCurveQuotient(D, N, W);
printf "=== X_0^%o(%o) / %o : genus %o, cut %o, %o points, %o digits  [oldforms included]\n",
       D, N, Sort(SetToSequence(W)), g, cutv, npts, PREC;
require g ge 3 : "Noether's count needs genus at least 3";
require IsSquarefree(D*N) and GCD(D, N) eq 1 : "the construction assumes D*N squarefree";
DN := D*N;
nmax := Max(60, Ceiling(DN*(cutv/(4*Pi(RR)))^2) + 20);

// ---------------------------------------------------------------- orders and Atkin-Lehner elements
A := QuaternionAlgebra(D);
// QuaternionOrder(A, N) depends on Magma's random state (a conjugate order each time the state
// differs), which moves the evaluation points on the curve.  Two runs of this script with the
// same preamble agree; SEED:=s makes the order reproducible against any other script that seeds
// at the same place.
if assigned SEED then SetSeed(StringToInteger(SEED)); end if;
ON := QuaternionOrder(A, N); O1 := MaximalOrder(ON);
mrep := FuchsianMatrixRepresentation(A : Precision := PREC);

function SuperOrder(ON, O1, M)
    A := Algebra(ON);
    return QuaternionOrder([ A ! x : x in Basis(ON) ] cat [ A ! (M*x) : x in Basis(O1) ]);
end function;

// a basis of J_m = {x in O_N : trd(x y) = 0 mod m for all y in O_N}
function TraceIdealBasis(ON, m)
    A := Algebra(ON);
    bas := [ A ! x : x in Basis(ON) ];
    T := Matrix(Integers(), 4, 4, [ Trace(bas[i]*bas[j]) : i, j in [1..4] ]);
    K := Kernel(Matrix(Integers(m), T));
    rows := [ [ Integers() ! c : c in Eltseq(K.i) ] : i in [1..Dimension(K)] ]
            cat [ Eltseq(m*IdentityMatrix(Integers(), 4)[i]) : i in [1..4] ];
    B := BasisMatrix(Lattice(Matrix(Integers(), rows)));
    assert Abs(Determinant(B)) eq m^2;
    return [ &+[ B[i][j]*bas[j] : j in [1..4] ] : i in [1..4] ];
end function;

// an element of O_N of reduced norm exactly m normalising O_N: a short vector of J_m for the
// Frobenius majorant of the real splitting
function ALElt(ON, m, mrep)
    A := Algebra(ON);
    jb := TraceIdealBasis(ON, m);
    R30 := RealField(30);
    Ms := [ ChangeRing(mrep(b), R30) : b in jb ];
    G := Matrix(Integers(), 4, 4, [ Round(10^6 * (&+[ Ms[a][r][c]*Ms[b][r][c] : r, c in [1..2] ])/2) : a, b in [1..4] ]);
    for a in [1..4] do for b in [a+1..4] do G[b][a] := G[a][b]; end for; G[a][a] +:= 1; end for;
    L := LatticeWithGram(G);
    bound := 10^6 * m * 4;
    while bound le 10^6 * m * 10^8 do
        for v in ShortVectors(L, bound) do
            x := &+[ (Integers() ! v[1][i]) * jb[i] : i in [1..4] ];
            if Norm(x) ne m then continue; end if;
            if forall{ b : b in Basis(ON) | x*(A!b)*x^-1 in ON } then return x; end if;
        end for;
        bound *:= 4;
    end while;
    error "no Atkin-Lehner element of norm", m;
end function;

function ALMul(u, v) gg := GCD(u, v); return (u div gg) * (v div gg); end function;
function ALGens(W)
    span := {Integers() | 1}; gens := [];
    for m in Sort(SetToSequence(W)) do
        if m in span then continue; end if;
        Append(~gens, m);
        span := { ALMul(s, m) : s in span } join span;
    end for;
    return gens;
end function;

// --------------------------------------------------------------------------- the forms, by level
// Sel: one entry per D-new newform kept, <M, cs, T, chis, rej>, with cs the complex coefficients
// a_g(n) for n <= nmax, T the primes of N/M, chis the characters giving W-fixed eigenvectors and
// rej the others (kept only with CHECKW, and only for oldclasses).
Sel := [* *]; newcnt := 0; oldcnt := 0;
for M in Divisors(N) do
    T := PrimeDivisors(N div M);
    NFs := Newforms(CuspForms(Gamma0(D*M), k));
    for orb in NFs do
        KF := BaseRing(Parent(orb[1]));
        cj := (Degree(KF) gt 1) select Conjugates(KF.1 : Precision := PREC+10) else [];
        for idx in [1..#orb] do
            f := orb[idx];
            if Level(Parent(f)) ne D*M then continue; end if;
            qq := qExpansion(f, nmax+1);
            if Degree(KF) eq 1 then
                cs := [ CC ! (Rationals() ! Coefficient(qq, n)) : n in [1..nmax] ];
            else
                cs := [ CC | Evaluate(Polynomial(Rationals(), Eltseq(Coefficient(qq, n))), cj[idx]) : n in [1..nmax] ];
            end if;
            chis := []; rej := [];
            for bits in CartesianPower([1, -1], #T) do
                chi := [ bits[i] : i in [1..#T] ];
                ok := true;
                for m in W do
                    if m eq 1 then continue; end if;
                    e := CC ! 1;
                    for p in PrimeDivisors(m) do
                        if D mod p eq 0 then e *:= cs[p];
                        elif M mod p eq 0 then e *:= -cs[p];
                        else e *:= chi[Index(T, p)]; end if;
                    end for;
                    if Abs(e - 1) gt 0.01 then ok := false; break; end if;
                end for;
                if ok then Append(~chis, chi); else Append(~rej, chi); end if;
            end for;
            if M eq N or nchk eq 0 then rej := []; end if;
            if #chis eq 0 and #rej eq 0 then continue; end if;
            Append(~Sel, <M, cs, T, chis, rej>);
            if M eq N then newcnt +:= #chis; else oldcnt +:= #chis; end if;
        end for;
    end for;
end for;
printf "W-fixed differentials: %o from newforms of level %o, %o from oldclasses; genus says %o\n",
       newcnt, DN, oldcnt, g;
for s in Sel do
    printf "   level %o (M=%o)  a_2..a_6 = %o  T = %o  chi fixed = %o%o\n", D*s[1], s[1],
           [ RR ! Real(s[2][n]) : n in [2..6] ], s[3], s[4],
           #s[5] gt 0 select Sprintf("  (carried, not fixed: %o)", s[5]) else "";
end for;
require newcnt + oldcnt eq g :
    "the W-fixed count including oldforms does not equal the genus; the eigenvalue bookkeeping is wrong";

// ------------------------------------------------------------- setups per level, columns, w_S
Wtab := WkTable(k, cutv + 2, RR);
levels := Sort(SetToSequence({ Integers() | s[1] : s in Sel }));
dat := AssociativeArray(); formsAt := AssociativeArray(); colsAt := AssociativeArray();
for M in levels do
    OM := (M eq N) select ON else SuperOrder(ON, O1, M);
    assert Discriminant(OM) eq D*M and forall{ x : x in Basis(ON) | A ! x in OM };
    dat[M] := NelsonSetupT(D, M, k, PREC, cutv + 2, Wtab : OO := OM, MREP := mrep);
    formsAt[M] := [ s[2] : s in Sel | s[1] eq M ];
    // columns: <form index within the level, chi, fixed?>; the fixed ones first
    cols := []; i := 0;
    for s in Sel do
        if s[1] ne M then continue; end if;
        i +:= 1;
        for chi in s[4] do Append(~cols, <i, chi, true>); end for;
    end for;
    i := 0;
    for s in Sel do
        if s[1] ne M then continue; end if;
        i +:= 1;
        for chi in s[5] do Append(~cols, <i, chi, false>); end for;
    end for;
    colsAt[M] := cols;
end for;
wp := AssociativeArray();
for p in PrimeDivisors(N) do
    w := ALElt(ON, p, mrep);
    assert Norm(w) eq p and (A ! (w*w))/p in ON;
    wp[p] := w;
    printf "w_%o = %o\n", p, w;
end for;

// Column labels in the order ColumnValues produces them: fixed columns of every level first (in
// level order), then the carried rejected ones.  The first g entries are the W-fixed basis.
labels := []; rejlabels := [];
for M in levels do
    for c in colsAt[M] do
        lab := Sprintf("level %o form %o chi %o", D*M, c[1], c[2]);
        if c[3] then Append(~labels, lab); else Append(~rejlabels, lab); end if;
    end for;
end for;
assert #labels eq g;

// predicted eigenvalue of w_m on a column
function PredEig(M, cs, T, chi, m)
    e := CC ! 1;
    for p in PrimeDivisors(m) do
        if D mod p eq 0 then e *:= cs[p];
        elif M mod p eq 0 then e *:= -cs[p];
        else e *:= chi[Index(T, p)]; end if;
    end for;
    return e;
end function;

// Values of every column at z (second point), z1 fixed.  Returns fixed columns, rejected columns,
// lattice-point count.
function ColumnValues(z1, z)
    fixed := [CC | ]; rejv := [CC | ]; nt := 0;
    for M in levels do
        T := PrimeDivisors(N div M);
        subs := Subsets({ p : p in T });
        vS := AssociativeArray();
        for S in subs do
            w := &*([A ! 1] cat [ A ! wp[p] : p in S ]);
            Mw := ChangeRing(mrep(w), CC);
            jj := Mw[2][1]*z + Mw[2][2];
            zz := (Mw[1][1]*z + Mw[1][2])/jj;
            vals, n1 := NelsonEvalMulti(dat[M], z1, zz, formsAt[M], RR ! cutv, CC);
            nt +:= n1;
            ph := ComplexConjugate(jj)/jj;
            vS[S] := [ ph * v : v in vals ];
        end for;
        for c in colsAt[M] do
            tot := CC ! 0;
            for S in subs do
                sg := &*([1] cat [ c[2][Index(T, p)] : p in S ]);
                tot +:= sg * vS[S][c[1]];
            end for;
            if c[3] then Append(~fixed, tot); else Append(~rejv, tot); end if;
        end for;
    end for;
    return fixed, rejv, nt;
end function;

// --------------------------------------------------------------------------------------- points
sp := Sqrt(RR ! ArithmeticVolume(FuchsianGroup(ON)));
z1 := CC ! (0.1103 + 1.0411*Ic);
zs := [ z1 ] cat [ CC | (sp*0.73*Cos(RR ! (1.1*j^2 + 0.3*j)))
                        + (0.55 + sp*0.42*(1 + Sin(RR ! (0.7*j^2 + 1.9*j))))*Ic
                   : j in [1..npts-1] ];
vals := []; rvals := [];
t0 := Cputime(); totpts := 0;
for j in [1..npts] do
    v, rv, nt := ColumnValues(z1, zs[j]);
    Append(~vals, v); Append(~rvals, rv);
    totpts +:= nt;
    if j eq 1 then printf "first point: %o s, %o lattice points\n", RR ! Cputime(t0), nt; end if;
end for;
printf "%o evaluations of %o columns: %o s, %o lattice points total\n",
       npts, g, RR ! Cputime(t0), totpts;

// ------------------------------------------------------------------- W-invariance of every column
if nchk gt 0 then
    printf "W-invariance check: V(w_m z) conj(j)/j against V(z); the entry is measured/predicted\n";
    gens := ALGens(W);
    worst := RR ! 0;
    for m in gens do
        wm := ALElt(ON, m, mrep);
        Mw := ChangeRing(mrep(wm), CC);
        for j in [1..Min(nchk, npts)] do
            z := zs[j]; jj := Mw[2][1]*z + Mw[2][2]; zz := (Mw[1][1]*z + Mw[1][2])/jj;
            v2, rv2 := ColumnValues(z1, zz);
            ph := ComplexConjugate(jj)/jj;
            // predicted eigenvalues, in column order
            preds := []; rpreds := [];
            for M in levels do
                i := 0; fs := [ s : s in Sel | s[1] eq M ];
                for c in colsAt[M] do
                    s := fs[c[1]];
                    e := PredEig(M, s[2], s[3], c[2], m);
                    if c[3] then Append(~preds, e); else Append(~rpreds, e); end if;
                end for;
            end for;
            rat := [ v2[c]*ph/vals[j][c]/preds[c] : c in [1..g] ];
            rrat := [ rv2[c]*ph/rvals[j][c]/rpreds[c] : c in [1..#rv2] ];
            printf "   w_%o, point %o: fixed columns %o\n", m, j, [ C8 ! r : r in rat ];
            if #rrat gt 0 then
                printf "      carried (predicted %o): %o\n", [ C8 ! e : e in rpreds ], [ C8 ! r : r in rrat ];
            end if;
            for r in rat cat rrat do worst := Max(worst, Abs(r - 1)); end for;
        end for;
    end for;
    printf "   largest deviation of measured/predicted from 1 over all columns and generators: %o\n", RR ! worst;
end if;

// ------------------------------------------------------------------------------- Noether's count
mons := [ [a,b] : b in [a..g], a in [1..g] ];
M := Matrix(CC, npts, #mons,
    &cat[ [ vals[j][m[1]]*vals[j][m[2]] : m in mons ] : j in [1..npts] ]);
for i in [1..npts] do
    nm := Sqrt(&+[ Abs(M[i][c])^2 : c in [1..#mons] ]);
    for c in [1..#mons] do M[i][c] := M[i][c]/nm; end for;
end for;
for c in [1..#mons] do
    nm := Sqrt(&+[ Abs(M[i][c])^2 : i in [1..npts] ]);
    for i in [1..npts] do M[i][c] := M[i][c]/nm; end for;
end for;
MR := Matrix(RR, 2*npts, 2*#mons, [ RR ! 0 : i in [1..4*npts*#mons] ]);
for i in [1..npts] do for j in [1..#mons] do
    a := Real(M[i][j]); b := Imaginary(M[i][j]);
    MR[i][j] := a; MR[i][#mons+j] := -b; MR[npts+i][j] := b; MR[npts+i][#mons+j] := a;
end for; end for;
Sig := SingularValueDecomposition(MR);
r := Min(Nrows(MR), Ncols(MR));
sv := Reverse([ Abs(Sig[i][i]) : i in [1..r] ]);
sv := [ sv[i] : i in [1..#sv by 2] ];
printf "singular values of the %o x %o monomial matrix:\n  %o\n", npts, #mons,
       [ RR ! x : x in sv ];
printf "consecutive ratios: %o\n", [ RR ! (sv[i+1]/sv[i]) : i in [1..#sv-1] ];

best := 1; bestr := RR ! 1;
for i in [1..#sv-1] do
    if sv[i] gt 0 and sv[i+1]/sv[i] gt bestr then bestr := sv[i+1]/sv[i]; best := i; end if;
end for;
mingap := RR ! 10^4;
if sv[1] gt 10^(-4) then
    nrel := 0; ok := true;
else
    nrel := best; ok := bestr ge mingap;
end if;
printf "largest ratio jump %o after %o value(s); clean gap required %o\n", RR!bestr, best, mingap;
if not ok then
    printf "=> NO CLEAN GAP: the spectrum is a ramp, the count does not resolve\n";
    printf "VERDICT: INCONCLUSIVE\n";
else
    printf "=> quadric relations: %o   (hyperelliptic would be %o, non-hyperelliptic %o)\n",
           nrel, (g-1)*(g-2) div 2, (g-2)*(g-3) div 2;
    printf "VERDICT: %o\n",
       nrel eq (g-1)*(g-2) div 2 select "HYPERELLIPTIC"
       else (nrel eq (g-2)*(g-3) div 2 select "NOT hyperelliptic" else "INCONCLUSIVE");
end if;
printf "smallest singular value %o ; old columns %o of %o\n", RR ! sv[1], oldcnt, g;
printf "STAGE3O DONE\n";
quit;
