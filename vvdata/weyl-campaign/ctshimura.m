// GAP A: prove  rho = [m_s ct(theta_s) - m_1 ct(theta_1)] / (m_s - m_1).
//
// Today's local identity collapses the RHS completely:
//   m_s ct_s - m_1 ct_1 = (c_2^G/2C) prod(p-1) prod c_p^ram,
//   m_s - m_1           = prod(p-1)/(2C)
//   =>  RHS = c_2^G * prod_{p|DN} c_p^ram.
// And the LHS factors as  rho = c_2^L * prod_{p|DN} c_p^{L,p}.
// At p | D the algebra is RAMIFIED, and ramification is LOCAL -- it does not
// see the archimedean place -- so c_p^{L,p} = c_p^ram, the same anisotropic
// plane as the definite case.  Hence the identity is EQUIVALENT to
//
//              c_2^L  =  c_2^G .
//
// This checks (i) the local forms of the Shimura lattice at odd p, (ii) the
// 2-parts of both lattices, and (iii) the identity itself.
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/fastcosets.m";
load "vvdata/weyl-campaign/etared.m";
load "vvdata/weyl-campaign/ctlat.m";
PREC := 120; CC := ComplexField(PREC);
frac := func< r | r - Floor(r) >;

// local structure of an arbitrary Gram: isotropic counts per prime + 2-part type
report := procedure(tag, G, QSIGN)
    dn := Denominator(G^(-1));
    Ld := RSpaceWithBasis(ChangeRing(dn*G^(-1), Integers()));
    Lp := RSpaceWithBasis(ScalarMatrix(3, dn));
    A, toA := Ld/Lp;
    mods := Moduli(A); keep := [r : r in [1..#mods] | mods[r] gt 1];
    ms := [mods[r] : r in keep]; k := #ms;
    gens := [ A.(keep[r]) : r in [1..k] ];
    wg := [ ChangeRing(g@@toA, Rationals()) : g in gens ];
    QG := [ frac(QSIGN*(wg[r]*G, wg[r])/(2*dn^2)) : r in [1..k] ];
    BG := [ [ frac(QSIGN*(wg[r]*G, wg[s])/dn^2) : s in [1..k] ] : r in [1..k] ];
    Qc := function(c)
        s := &+[ Rationals() | c[r]^2*QG[r] : r in [1..k] ];
        for r in [1..k] do for s2 in [r+1..k] do s +:= c[r]*c[s2]*BG[r][s2]; end for; end for;
        return frac(s);
    end function;
    printf "  %-22o |A| = %-9o invariants %o\n", tag, &*ms, ms;
    for p in PrimeDivisors(&*ms) do
        JP := [ [Integers()|0 : r in [1..k]] ];
        for r in [1..k] do
            ar := Valuation(ms[r], p); step := ms[r] div p^ar; nxt := [];
            for c in JP do for t in [0..p^ar-1] do
                cc := c; cc[r] := (t*step) mod ms[r]; Append(~nxt, cc); end for; end for;
            JP := nxt;
        end for;
        if p eq 2 then
            nzs := [ j : j in JP | not forall{r : r in [1..k] | j[r] mod ms[r] eq 0} ];
            printf "      p=2: |J_2|=%o, 4*Qbar on the nonzero element(s) = %o\n",
                #JP, [ Integers()!(4*Qc(j)) mod 4 : j in nzs ];
        else
            nz := #[ j : j in JP | Qc(j) eq 0 ];
            printf "      p=%2o: |J_p|=%4o, %o isotropic -> %o\n", p, #JP, nz,
                (nz eq 1) select "ANISOTROPIC" else ((nz eq 2*p-1) select "HYPERBOLIC" else "?");
        end if;
    end for;
end procedure;

for DN in [1155, 330] do
    printf "\n===== DN = %o =====\n", DN;
    Ld := ShimuraCurveLattice(DN, 1);
    printf "SHIMURA CURVE LATTICE (indefinite, sig (1,2)), QSIGN = -1 as eis32s uses:\n";
    report(Sprintf("L(%o,1)", DN), ChangeRing(Ld`Q, Rationals()), -1);
    allp := PrimeDivisors(DN);
    s := allp[#allp]; Dp := DN div s;
    printf "GROSS LATTICES of the support (%o,1)+(%o,%o), QSIGN = -1:\n", Dp, Dp, s;
    report(Sprintf("G(%o,1)", Dp), grossGram(Dp, 1), -1);
    report(Sprintf("G(%o,%o)", Dp, s), grossGram(Dp, s), -1);
end for;
printf "\nDONE\n";
quit;
