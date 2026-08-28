// X_0^1155(1), Q4: the five preregistered candidate supports, by CONSTANT
// TERMS.  Same machinery as ctgate330.m, which reproduces the banked 330_1
// answer -2*theta(165,1) + 3*theta(165,2) to 115 digits with 0 of 3456 cosets
// missing.  Here the cuspidal basis never appears, so the level-4620 wall that
// stopped Magma for 11 h is simply absent.
//
// prereg-1155.md predicts support {(385,1),(385,3)} with weights (-1,+2) and
// the OTHER FOUR FAILING.  All five are tested: "several supports carry E_eis"
// is the informative outcome and cannot be seen by testing only the winner.
//
// allgross_330_1.log says
//     E_eis = -2 * GROSS(165,1) + 3 * GROSS(165,2)   (+ cusp terms)
// Under the constant-term formulation the cusp terms have NO constant term at
// any cusp, so they must drop out identically and the relation
//     ct(E_eis) = -2 * ct(theta(165,1)) + 3 * ct(theta(165,2))
// must hold EXACTLY, at all 1728 cosets.  That is a sharp check on the
// derivation, on ctlat.m, and on the sign convention QSIGN -- so both signs
// are tried and the gate reports which (if either) reproduces the banked
// weights.
//   magma -b vvdata/weyl-campaign/ctgate330.m
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/fastcosets.m";
load "vvdata/weyl-campaign/etared.m";
load "vvdata/weyl-campaign/ctlat.m";

PREC := 120;
CC := ComplexField(PREC); ii := CC.1; pi := Pi(CC);
ee := func< z | Exp(2*pi*ii*z) >;

// TWO DIFFERENT LEVELS, and they must not be conflated (the handoff records
// this error costing a preregistered prediction):
//   * the ETA level Meta = 2*DN = 660 indexes the eta monomials (dsv);
//   * the weight-3/2 FORM level lev = 4*DN = 1320 indexes the COSETS.
// At 1155 they coincide (odd DN => Meta = 4DN), which is exactly why this only
// shows up in the gate.
DN := 1155; lev := 4*DN;                   // 4620: the form level
Meta := IsOdd(DN) select 4*DN else 2*DN;   // 4620: odd DN, so they coincide
ds := Divisors(Meta); nd := #ds;
words := [ VVSTWord(g) : g in fastCosetReps(lev) ];
nw := #words;
printf "1155_1 Q4: form level %o (%o cosets), eta level %o (%o divisors)\n",
    lev, nw, Meta, nd;

// ---- E_eis at 330_1, exactly the data allgross330_1.m uses -----------------
src := Read("vvdata/weyl-campaign/allgross1155_1.m");
// pull dsv, rs, cs out of the generated driver so this cannot drift from it
getblk := function(tag)
    i := Index(src, tag);
    error if i eq 0, "tag not found", tag;
    j := i;
    depth := 0; started := false;
    while true do
        ch := src[j];
        if ch eq "[" then depth +:= 1; started := true; end if;
        if ch eq "]" then depth -:= 1; end if;
        if started and depth eq 0 then break; end if;
        j +:= 1;
    end while;
    return src[i..j];
end function;
dsv := eval(getblk("dsv := ")[8..#getblk("dsv := ")]);
error if dsv ne ds, "divisor list mismatch: driver dsv is not Divisors(eta level)";
rsblk := getblk("rs := "); rs := eval(rsblk[7..#rsblk]);
csblk := getblk("cs := "); cs := eval(csblk[7..#csblk]);
printf "E_eis: %o monomials\n", #rs;

// ---- a0at: constant term of an eta monomial at a coset (as in eis32s) ------
triang := function(g, d)
    g2 := Matrix(Integers(), 2, 2, [d*g[1][1], d*g[1][2], g[2][1], g[2][2]]);
    c1 := g2[1][1]; c2 := g2[2][1]; h := GCD(c1, c2);
    p1 := c1 div h; p2 := c2 div h; gg, u, v := XGCD(p1, p2);
    gd := Matrix(Integers(), 2, 2, [p1, -v, p2, u]); sd := gd^(-1) * g2;
    a := sd[1][1]; b := sd[1][2]; e := sd[2][2];
    if a lt 0 then a := -a; b := -b; e := -e; end if;
    if e lt 0 then sd := -sd; a := sd[1][1]; b := sd[1][2]; e := sd[2][2]; end if;
    return a, b, e;
end function;
slashdata := function(word, tau)
    z := tau; factor := CC!1;
    for i := #word to 1 by -1 do
        if word[i][1] eq "S" then factor /:= Sqrt(z); z := -1/z; else z := z + word[i][2]; end if;
    end for;
    return factor, z;
end function;
tau0 := CC!0.31 + CC!1.31*ii; tau1 := CC!(-0.57) + CC!1.73*ii;
SS := PowerSeriesRing(CC); t := SS.1;

ctE := [ CC | 0 : wi in [1..nw] ];
for wi->w in words do
    g := VVWordMatrix(w);
    tri := [ ]; for d in ds do a,b,e := triang(g,d); Append(~tri, <a,b,e>); end for;
    W := LCM([ tri[i][3] : i in [1..nd] ]);
    fac0, z0 := slashdata(w, tau0); fac1, z1 := slashdata(w, tau1);
    rat0 := [ CC | ]; rat1 := [ CC | ];
    for i->d in ds do
        a := tri[i][1]; b := tri[i][2]; e := tri[i][3];
        s0 := etaRed((a*tau0 + b)/e) * ee(-(a*tau0 + b)/(24*e));
        s1 := etaRed((a*tau1 + b)/e) * ee(-(a*tau1 + b)/(24*e));
        Append(~rat0, etaRed(d*z0)/s0); Append(~rat1, etaRed(d*z1)/s1);
    end for;
    tot := CC!0;
    for ei in [1..#rs] do
        rE := rs[ei];
        L := &+[ Integers() | rE[i]*tri[i][1]*(W div tri[i][3]) : i in [1..nd] ];
        if L gt 0 then continue; end if;
        error if L ne 0, "unexpected depth > 1 at 1155", wi, ei, L;
        p0 := CC!1; p1 := CC!1;
        for i in [1..nd] do
            if rE[i] eq 0 then continue; end if;
            p0 *:= rat0[i]^(rE[i]); p1 *:= rat1[i]^(rE[i]);
        end for;
        k0 := fac0^3 * p0; k1 := fac1^3 * p1;
        error if Abs(k0 - k1) gt 10^(-25)*Maximum(Abs(k0), 1), "kappa not constant", wi, ei;
        tot +:= (CC!cs[ei]) * k0;      // c0 = 1 since L = 0
    end for;
    ctE[wi] := tot;
end for;
printf "ct(E_eis) built\n";

// ---- Q4: each candidate support, QSIGN = -1 as pinned by the 330_1 gate ----
QSIGN := -1;
cands := [
    <"385", 385, 1, 3,   -1, 2>,      // predicted winner, weights (-1,+2)
    <"105", 105, 1, 11,  0, 0>,
    <"11",  11,  1, 105, 0, 0>,
    <"231", 231, 1, 5,   0, 0>,
    <"165", 165, 1, 7,   0, 0>
];
printf "\n==== Q4 at 1155: ct(E_eis) in span{ct(theta1), ct(theta2)} ? ====\n";
for c in cands do
    Dp := c[2]; R1 := c[3]; R2 := c[4];
    G1 := grossGram(Dp, R1); G2 := grossGram(Dp, R2);
    c1v, n1 := ctTheta(G1, words, CC, QSIGN);
    c2v, n2 := ctTheta(G2, words, CC, QSIGN);
    a11 := &+[ c1v[i]*ComplexConjugate(c1v[i]) : i in [1..nw] ];
    a12 := &+[ c1v[i]*ComplexConjugate(c2v[i]) : i in [1..nw] ];
    a22 := &+[ c2v[i]*ComplexConjugate(c2v[i]) : i in [1..nw] ];
    b1  := &+[ ctE[i]*ComplexConjugate(c1v[i]) : i in [1..nw] ];
    b2  := &+[ ctE[i]*ComplexConjugate(c2v[i]) : i in [1..nw] ];
    det := a11*a22 - a12*ComplexConjugate(a12);
    printf "\nSUPPORT (%o,%o)+(%o,%o):  |A| = %o, %o\n", Dp, R1, Dp, R2, n1, n2;
    if Abs(det) lt 10^(-60) then printf "  DEGENERATE (theta columns dependent)\n"; continue; end if;
    x1 := (b1*a22 - b2*a12)/det;
    x2 := (a11*b2 - ComplexConjugate(a12)*b1)/det;
    res := Sqrt(&+[ Abs(ctE[i] - x1*c1v[i] - x2*c2v[i])^2 : i in [1..nw] ]);
    nrm := Sqrt(&+[ Abs(ctE[i])^2 : i in [1..nw] ]);
    inspan := res lt 10^(-20)*nrm;
    printf "  w1 = %o\n  w2 = %o\n", ComplexField(24)!x1, ComplexField(24)!x2;
    printf "  residual = %o  relative = %o\n", RealField(10)!res, RealField(10)!(res/nrm);
    printf "  E_eis IN SPAN: %o\n", inspan;
    if inspan then
        printf "  weight sum = %o   |w2/w1| = %o\n",
            ComplexField(20)!(x1 + x2), RealField(20)!(Abs(x2)/Abs(x1));
    end if;
end for;

// ---- IS THE TEST DEGENERATE?  ----------------------------------------------
// Four supports passing with exactly their predicted weights is either a real
// finding or an artifact of the ct(theta) vectors spanning a tiny space.  If
// every pair spans the SAME plane, the test cannot distinguish supports at all
// and "four supports carry E_eis" says nothing.  Measure the rank.
printf "\n==== degeneracy check: rank of the ct(theta) system ====\n";
allp := [ <385,1>,<385,3>,<105,1>,<105,11>,<11,1>,<11,105>,<231,1>,<231,5>,<165,1>,<165,7> ];
vecs := [ ]; nms := [ ];
for pr in allp do
    v := ctTheta(grossGram(pr[1], pr[2]), words, CC, QSIGN);
    Append(~vecs, v); Append(~nms, Sprintf("(%o,%o)", pr[1], pr[2]));
end for;
// numerical rank by Gram-matrix eigenvalues (vectors are complex, length nw)
np := #vecs;
Gm := Matrix(CC, np, np, [ &+[ vecs[i][t]*ComplexConjugate(vecs[j][t]) : t in [1..nw] ]
                           : j in [1..np], i in [1..np] ]);
// normalise each vector so the scale is comparable
nrms := [ Sqrt(Re(Gm[i][i])) : i in [1..np] ];
Gn := Matrix(CC, np, np, [ Gm[i][j]/(nrms[i]*nrms[j]) : j in [1..np], i in [1..np] ]);
printf "pairwise |<vi,vj>| (normalised), 1.000 means PARALLEL:\n      ";
for j in [1..np] do printf "%9o", nms[j]; end for; printf "\n";
for i in [1..np] do
    printf "%9o", nms[i];
    for j in [1..np] do printf "  %7o", RealField(5)!Abs(Gn[i][j]); end for;
    printf "\n";
end for;
// rank via eigenvalues of the normalised Gram
ev := Eigenvalues(Matrix(RealField(60), np, np, [ Re(Gn[i][j]) : j in [1..np], i in [1..np] ]));
evs := Sort([ RealField(12)!e[1] : e in ev ]);
printf "eigenvalues of the normalised Gram: %o\n", evs;
printf "numerical rank (eig > 1e-30): %o of %o\n",
    #[ e : e in evs | e gt 10^(-30) ], np;
printf "DONE\n";
quit;
