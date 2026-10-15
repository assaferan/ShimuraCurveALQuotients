// Does the Hauptmodul of a star curve computed from the expansions agree with the published one?
//
// The star curve X_0^D(N)* is the quotient by the FULL Atkin-Lehner group, so this exercises the
// Atkin-Lehner layer as well as the base expansions.  The comparison is coordinate-free: a computed
// Hauptmodul is defined only up to a Mobius change, so three CM points fix the coordinate and every
// further CM value is then a prediction, tested as a cross-ratio against the published one.
//
// At D = 10 the published values are Elkies, "Shimura curve computations" (arXiv:math/0005160),
// Sec. 4 and Table 3, as transcribed with its corrections in tests/CMFieldsElkies.m on main; the
// four elliptic points sit at t = 0, oo, 2, 27 for discriminants -3, -8, -20, -40.  These are
// INDEPENDENT of this repository's own model files.
//
//   magma -b DISC:=10 LEVEL:=1 PR:=140 TGT:=16 WT:=12 BOX:=3 EBOX:=8 starcmp.m < /dev/null
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vw3.m";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vwal.m";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/cmpts.m";

D := StringToInteger(DISC); N := StringToInteger(LEVEL);
PREC := StringToInteger(PR); TARGET := StringToInteger(TGT);
k := StringToInteger(WT); BOXSZ := StringToInteger(BOX); EB := StringToInteger(EBOX);
DN := D*N; CC := ComplexField(PREC);

// ---- published values: discriminant -> Hauptmodul value of the STAR curve (Infinity() allowed).
// Elkies, "Shimura curve computations" (arXiv:math/0005160): Sec. 3.1 and Table 1 for D = 6,
// Sec. 4 and Table 3 for D = 10, Sec. 5.1 and Table 5 for D = 14, Sec. 5.2 and Table 6 for D = 15,
// transcribed with their corrections in tests/CMFieldsElkies.m on main.  Only FUNDAMENTAL
// discriminants are listed: a non-fundamental one needs the embedding checked for optimality,
// since an element generating the smaller order fixes the point of the larger one.
QQ := Rationals();
ELK := AssociativeArray();
ELK[6]  := [* <-3, Infinity()>, <-4, QQ!1>, <-24, QQ!0>, <-84, -169/27>, <-40, 2312/125>,
              <-51, -1377/1024>, <-19, 3211/1024>, <-120, 5776/3375>, <-52, 6877/15625>,
              <-132, 13689/15625>, <-168, -701784/15625>, <-43, 21250987/16000000>,
              <-88, 15545888/20796875>, <-312, 27008742384/27680640625> *];
ELK[10] := [* <-3, QQ!0>, <-8, Infinity()>, <-20, QQ!2>, <-40, QQ!27>, <-52, -54/25>,
              <-120, -27/49>, <-35, 64/7>, <-43, 1728/1225>, <-88, 3375/98>,
              <-115, 13824/3887>, <-280, 35937/7406>, <-67, -216000/8281> *];
ELK[14] := [* <-4, Infinity()>, <-8, QQ!0>, <-11, QQ!-1>, <-35, QQ!-5>, <-84, 4/9>, <-91, -13/81>,
              <-43, -25/81>, <-168, 200/9>, <-88, -200/81>, <-67, -1225/81>, <-280, -845/1296>,
              <-148, 1225/5184> *];
ELK[15] := [* <-15, QQ!81>, <-3, Infinity()>, <-7, QQ!-27>, <-40, 27/2>, <-43, -27/16>,
              <-195, 81/16>, <-120, -243/2>, <-115, 621/16>, <-123, 2025/16>, <-67, -3267/16> *];
error if not IsDefined(ELK, D), Sprintf("no published table recorded for D = %o", D);
elk := ELK[D];
// the first three FINITE published values fix the coordinate; everything after is a prediction
refs := [ e[1] : e in elk | not (e[2] cmpeq Infinity()) ][1..3];

t0 := Cputime();
A := QuaternionAlgebra(D); O := QuaternionOrder(A, N);
G := FuchsianGroup(O);
H := UpperHalfPlane();
_, DD, FD := RecentredDisc(G, H, PREC, 3);
rho := Max([Abs(ComplexValue(v)) : v in FD]);
p := ComplexValue(DiscToPlane(H, DD ! 0) : Precision := PREC);
mrep := FuchsianMatrixRepresentation(A : Precision := PREC);
NT := TermsFor(rho, TARGET);
tol := ZeroLevel(rho, NT : prec := PREC);
Q := 3*NT;
printf "references: d = %o with published values %o\n", refs, [ [e[2] : e in elk | e[1] eq d][1] : d in refs ];
printf "D=%o N=%o vol %o rho %o sides %o | NT %o tol %o prec %o (setup %o s)\n",
    D, N, RealField(5)!ArithmeticVolume(G), RealField(6)!rho, #FD, NT,
    RealField(3)!tol, PREC, RealField(4)!Cputime(t0);

// ---- samples: Gamma, plus the Atkin-Lehner generators of the full group ----
Wfull := {Integers() | d : d in Divisors(DN) | GCD(d, DN div d) eq 1};
gens := ALGens(Wfull, DN);
mus := [];
for m in gens do
    ok, x := ALElement(O, m, BOXSZ);
    error if not ok, Sprintf("no Atkin-Lehner element of norm %o in a box of %o", m, BOXSZ);
    Append(~mus, x);
end for;
t1 := Cputime();
samp := ReducedSamples(G, O, DD, mrep, p, rho, Q, CC);
nG := #samp;
alsamp, nbad := ALSamples(G, O, DD, mrep, p, rho, Q, CC, mus);
samp cat:= alsamp;
printf "samples: %o from Gamma, %o from Atkin-Lehner (%o dropped), total %o (%o s)\n",
    nG, #alsamp, nbad, #samp, RealField(4)!Cputime(t1);

// ---- the weight-k space: expect dimension 2, and the ratio of a basis is a Hauptmodul ----
X := CreateShimuraQuot(D, N, Wfull);
g, es := QuotientSignature(X);
expdim := DimMkSig(k, g, es);
printf "star signature (%o; %o), expected dim M_%o = %o\n", g, es, k, expdim;

t2 := Cputime();
M := AutSystem(samp, k, NT, CC);
sc := [ Sqrt(&+[ Abs(M[i][j])^2 : i in [1..Nrows(M)] ]) : j in [1..NT] ];
for j in [1..NT] do for i in [1..Nrows(M)] do M[i][j] := M[i][j]/sc[j]; end for; end for;
MR := RealForm(M);
Sig, L, Rt := SingularValueDecomposition(MR);
r := Min(Nrows(MR), Ncols(MR));
svall := Reverse([ Abs(Sig[i][i]) : i in [1..r] ]);
sv := [ svall[i] : i in [1..#svall by 2] ];
gotdim, gap, glo, ghi := GapDimRel(sv, 1e-4);
printf "dim M_%o: expected %o got %o%o | gap %o between %o and %o | smallest %o (%o s)\n",
    k, expdim, gotdim, gotdim eq expdim select "" else "  MISMATCH", RealField(3)!gap,
    RealField(3)!glo, RealField(3)!ghi, [RealField(3)!s : s in sv[1..Min(6,#sv)]],
    RealField(4)!Cputime(t2);
error if gotdim ne expdim, "weight-k dimension is wrong; the Hauptmodul would be meaningless";

function CoeffVec(row) return [ CC | (Rt[row][i] + CC.1*Rt[row][NT+i])/sc[i] : i in [1..NT] ]; end function;
b1 := CoeffVec(r); b2 := CoeffVec(r-2);          // complex singular values appear twice
F := func< bs, w | &+[ bs[n+1]*w^n : n in [0..NT-1] ] >;
Haupt := func< w | F(b1, w)/F(b2, w) >;

// ---- CM points, and the cross-ratio test ----
vals := AssociativeArray(); info := AssociativeArray();
for e in elk do
    d := e[1];
    ok, w, mu, tn := CMDiscPoint(O, d, DD, mrep, p, CC, G);
    if not ok then printf "  d = %-5o no optimal embedding found\n", d; continue; end if;
    vals[d] := Haupt(w); info[d] := <Abs(w), tn>;
end for;
printf "CM points located: %o of %o\n", #Keys(vals), #elk;

// A computed Hauptmodul is defined only up to a Mobius change, so the test is whether ONE Mobius
// map carries the published values onto the computed ones at every CM point at once.  Writing
// s = (a t + b)/(c t + d) gives one linear equation a*t + b - c*t*s - d*s = 0 per point (and
// a - c*s = 0 when the published value is infinite), so n points give n equations in 4 unknowns
// defined up to scale: three points fix the coordinate and every further point is a prediction.
// The smallest singular value is the residual of the fit, the next one its conditioning.
ds := Sort([ d : d in Keys(vals) ]);
rows := [];
for d in ds do
    tv := [e[2] : e in elk | e[1] eq d][1];
    sv1 := vals[d];
    if tv cmpeq Infinity() then
        r := [ CC | 1, 0, -sv1, 0 ];
    else
        tc := CC ! (QQ ! tv);
        r := [ CC | tc, 1, -tc*sv1, -sv1 ];
    end if;
    nm := Sqrt(&+[ Abs(x)^2 : x in r ]);
    Append(~rows, [ x/nm : x in r ]);
end for;
MM := Matrix(CC, #rows, 4, &cat rows);
svM := SingVals(MM);
printf "\nMobius fit over %o CM points: singular values %o\n", #rows, [RealField(4)!x : x in svM];
printf "  residual %o, next %o, ratio %o\n", RealField(4)!svM[1], RealField(4)!svM[2],
    RealField(4)!(svM[2]/svM[1]);
vec := SmallestVec(MM);
a := vec[1]; b := vec[2]; c := vec[3]; dd := vec[4];
printf "\n%-7o %-22o %-24o %o\n", "d", "published t", "residual |s - M(t)|", "relative";
nok := 0;
for d in ds do
    tv := [e[2] : e in elk | e[1] eq d][1];
    if tv cmpeq Infinity() then pred := a/c; else tc := CC ! (QQ ! tv); pred := (a*tc + b)/(c*tc + dd); end if;
    res := Abs(vals[d] - pred);
    rel := Abs(vals[d]) gt 0 select res/Abs(vals[d]) else res;
    if rel lt 1e-8 then nok +:= 1; end if;
    printf "%-7o %-22o %-24o %o\n", d, tv cmpeq Infinity() select "oo" else Sprint(tv),
        RealField(4)!res, RealField(4)!rel;
end for;
printf "\nSTARCMP D=%o N=%o weight %o: %o of %o published CM values fitted by one Mobius map\n", D, N, k, nok, #ds;
quit;
