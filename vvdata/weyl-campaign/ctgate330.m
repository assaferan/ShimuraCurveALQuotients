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
DN := 330; lev := 4*DN;                    // 1320: the form level
Meta := IsOdd(DN) select 4*DN else 2*DN;   // 660: the eta level
ds := Divisors(Meta); nd := #ds;
words := [ VVSTWord(g) : g in fastCosetReps(lev) ];
nw := #words;
printf "330_1 gate: form level %o (%o cosets), eta level %o (%o divisors)\n",
    lev, nw, Meta, nd;

// ---- E_eis at 330_1, exactly the data allgross330_1.m uses -----------------
src := Read("vvdata/weyl-campaign/allgross330_1.m");
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
        error if L ne 0, "unexpected depth > 1 at 330_1", wi, ei, L;
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

// ---- the two Gross thetas of the banked support, both sign conventions -----
for QSIGN in [1, -1] do
    printf "\n---- QSIGN = %o ----\n", QSIGN;
    G1 := grossGram(165, 1); G2 := grossGram(165, 2);
    c1v, n1 := ctTheta(G1, words, CC, QSIGN);
    c2v, n2 := ctTheta(G2, words, CC, QSIGN);
    printf "  |A| = %o and %o\n", n1, n2;
    // least squares for ct(E) = x1 c1v + x2 c2v
    a11 := &+[ c1v[i]*ComplexConjugate(c1v[i]) : i in [1..nw] ];
    a12 := &+[ c1v[i]*ComplexConjugate(c2v[i]) : i in [1..nw] ];
    a22 := &+[ c2v[i]*ComplexConjugate(c2v[i]) : i in [1..nw] ];
    b1  := &+[ ctE[i]*ComplexConjugate(c1v[i]) : i in [1..nw] ];
    b2  := &+[ ctE[i]*ComplexConjugate(c2v[i]) : i in [1..nw] ];
    det := a11*a22 - a12*ComplexConjugate(a12);
    if Abs(det) lt 10^(-60) then printf "  degenerate\n"; continue; end if;
    x1 := (b1*a22 - b2*a12)/det;
    x2 := (a11*b2 - ComplexConjugate(a12)*b1)/det;
    res := Sqrt(&+[ Abs(ctE[i] - x1*c1v[i] - x2*c2v[i])^2 : i in [1..nw] ]);
    nrm := Sqrt(&+[ Abs(ctE[i])^2 : i in [1..nw] ]);
    printf "  x1 = %o\n  x2 = %o\n", ComplexField(20)!x1, ComplexField(20)!x2;
    printf "  residual = %o   |ct(E)| = %o   relative = %o\n",
        RealField(10)!res, RealField(10)!nrm, RealField(10)!(res/nrm);
    printf "  BANKED: -2 and +3\n";
    printf "  GATE: %o\n",
        (Abs(x1 + 2) lt 10^(-20) and Abs(x2 - 3) lt 10^(-20) and res lt 10^(-20)*nrm)
        select "PASS" else "no";
    if QSIGN eq -1 then
        // localize: use the BANKED weights exactly and see which cosets miss
        printf "  --- per-coset residual with the exact banked weights -2, +3 ---\n";
        bad := [ ]; worst := RealField(20)!0;
        for i in [1..nw] do
            d := Abs(ctE[i] - (CC!-2)*c1v[i] - (CC!3)*c2v[i]);
            if d gt 10^(-20) then Append(~bad, <i, d>); end if;
            if d gt worst then worst := RealField(20)!d; end if;
        end for;
        printf "  %o of %o cosets miss; worst %o\n", #bad, nw, worst;
        for j in [1..Minimum(#bad, 12)] do
            wi := bad[j][1]; g := VVWordMatrix(words[wi]); gi := g^(-1);
            printf "    coset %o: |miss| = %o  c = %o  gcd(c,lev) = %o  a(inv) = %o  #S = %o\n",
                wi, RealField(8)!bad[j][2], g[2][1], GCD(g[2][1] mod lev, lev), gi[1][1],
                #[ t : t in words[wi] | t[1] eq "S" ];
        end for;
        // how do the misses distribute by cusp class?
        cls := AssociativeArray();
        for b in bad do
            g := VVWordMatrix(words[b[1]]); c := GCD(g[2][1] mod lev, lev);
            if IsDefined(cls, c) then cls[c] +:= 1; else cls[c] := 1; end if;
        end for;
        printf "  misses by cusp class gcd(c,lev): %o\n",
            [ <c, cls[c]> : c in Sort([ x : x in Keys(cls) ]) ];
        // and the total count per class for comparison
        tot := AssociativeArray();
        for i in [1..nw] do
            g := VVWordMatrix(words[i]); c := GCD(g[2][1] mod lev, lev);
            if IsDefined(tot, c) then tot[c] +:= 1; else tot[c] := 1; end if;
        end for;
        printf "  total per class:                %o\n",
            [ <c, tot[c]> : c in Sort([ x : x in Keys(cls) ]) ];
    end if;
end for;
printf "DONE\n";
quit;
