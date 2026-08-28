// GATE for the constant-term route, at 330_1 where the answer is banked.
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
DN := 210; lev := 4*DN;                    // 840: the form level
Meta := IsOdd(DN) select 4*DN else 2*DN;   // 420: the eta level
ds := Divisors(Meta); nd := #ds;
words := [ VVSTWord(g) : g in fastCosetReps(lev) ];
nw := #words;
printf "210_1 control: form level %o (%o cosets), eta level %o (%o divisors)\n",
    lev, nw, Meta, nd;

// ---- E_eis at 330_1, exactly the data allgross330_1.m uses -----------------
src := Read("vvdata/weyl-campaign/allgross210_1.m");
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
        error if L ne 0, "unexpected depth > 1 at 210_1", wi, ei, L;
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

// ---- THE CONTROL: does the constant-term test DISCRIMINATE at 330_1? -------
// At 1155 four of five candidate supports passed.  Before reading that as
// "the support is a convention", check the method on a base where the banked
// q-expansion run says the competitor IS excluded (330_1: rank 135 vs 136).
// The analogous family here is {(D',1), (D', 330/D')} over the definite D'
// (odd number of primes).  Banked winner: D' = 165 with (-2, +3).
//   * if ONLY 165 passes -> the test discriminates, and the 1155 four-way pass
//     is a genuine property of that base;
//   * if SEVERAL pass -> the constant-term test is weaker than the q-expansion
//     test, and the 1155 result cannot be read as a statement about supports.
QSIGN := -1;
printf "\n==== CONTROL at 210_1: all eight candidate supports ====\n";
npass := 0;
for Dp in [2, 3, 5, 7, 30, 42, 70, 105] do
    R2 := 210 div Dp;
    G1 := grossGram(Dp, 1); G2 := grossGram(Dp, R2);
    c1v, n1 := ctTheta(G1, words, CC, QSIGN);
    c2v, n2 := ctTheta(G2, words, CC, QSIGN);
    a11 := &+[ c1v[i]*ComplexConjugate(c1v[i]) : i in [1..nw] ];
    a12 := &+[ c1v[i]*ComplexConjugate(c2v[i]) : i in [1..nw] ];
    a22 := &+[ c2v[i]*ComplexConjugate(c2v[i]) : i in [1..nw] ];
    b1  := &+[ ctE[i]*ComplexConjugate(c1v[i]) : i in [1..nw] ];
    b2  := &+[ ctE[i]*ComplexConjugate(c2v[i]) : i in [1..nw] ];
    det := a11*a22 - a12*ComplexConjugate(a12);
    if Abs(det) lt 10^(-60) then
        printf "SUPPORT (%o,1)+(%o,%o): DEGENERATE\n", Dp, Dp, R2; continue;
    end if;
    x1 := (b1*a22 - b2*a12)/det;
    x2 := (a11*b2 - ComplexConjugate(a12)*b1)/det;
    res := Sqrt(&+[ Abs(ctE[i] - x1*c1v[i] - x2*c2v[i])^2 : i in [1..nw] ]);
    nrm := Sqrt(&+[ Abs(ctE[i])^2 : i in [1..nw] ]);
    inspan := res lt 10^(-20)*nrm;
    if inspan then npass +:= 1; end if;
    printf "SUPPORT (%o,1)+(%o,%o): rel resid %o  IN SPAN: %o", Dp, Dp, R2,
        RealField(8)!(res/nrm), inspan;
    if inspan then
        printf "   w = (%o, %o)  sum %o  ratio %o",
            ComplexField(12)!x1, ComplexField(12)!x2,
            ComplexField(12)!(x1+x2), RealField(12)!(Abs(x2)/Abs(x1));
    end if;
    printf "\n";
end for;
printf "\nCONTROL VERDICT: %o of 8 supports pass at 210_1 (banked tested only (7,1)+(7,30), excluded)\n", npass;
printf "  => the constant-term test %o\n",
    npass eq 1 select "DISCRIMINATES" else "does NOT discriminate";
printf "DONE\n";
quit;
