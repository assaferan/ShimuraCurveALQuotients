// Isotropic-sum via overlattices (kappa = 2, even).
//
// For N_beta = 2 the only character mod 2 is trivial (even => contributes) but imprimitive, so Thm 1.4
// needs Prop 1.3's oldform lifting. With psi trivial (N_psi = 1): N_0 = 1, N_0' = N_beta = 2, 2*beta = 0,
//     E_{A,beta} = E_{B_1,0}^ - E_{A,0},   B_1 = disc form of the even OVERLATTICE L_beta = L + Z*beta.
// Hence for a base with two nonzero isotropic beta (15_2, 21_2):
//     G = sum_{beta iso} E_{A,beta} = E_{B_1,0}^ + E_{B_2,0}^ - E_{A,0}.
// At gamma = 0 the lift is trivial ((f^)_0 = f_{0+H}), so we just need each overlattice's beta = 0
// Eisenstein coefficient at its zero class -- the SAME density formula, different lattice.
//
// Per-lattice coefficient: b^X_0(n) = (96 sqrt2 / sqrt|det X|) * cond * (h/w) * (en/ed)_X * g_X(n).
// Only sqrt|det X|, S_c and g_X depend on the lattice (cond, h, w depend on n and D_0 only, and D_0 is
// unchanged since det scales by a square: 1800 -> 450 keeps D = -n*(square)).
// So relative to L the overlattice terms carry a factor sqrt(|det L| / |det X|) = 2 (index 2).
//
// TEST: this must PRESERVE the Table-45-validated 15_2 form[-1] -> 4. If it breaks it, the isotropic-sum
// route is dead for a second, independent reason.

AttachSpec("ShimuraQuotients.spec");
fh := Open("over_out.txt", "w");
emit := procedure(s) fprintf fh, "%o\n", s; Flush(fh); end procedure;

rep_count := function(G, gam, n, p, w)
    pw := p^w;
    if w le 0 then return 1, pw; end if;
    Gr := ChangeRing(G, Rationals());
    Gg := gam*Gr; cst := (gam*Gr, gam)/2 + n;
    assert IsIntegral(cst);
    for i in [1..3] do assert IsIntegral(Gg[i]); end for;
    h11 := Integers()!(G[1][1]/2); h22 := Integers()!(G[2][2]/2); h33 := Integers()!(G[3][3]/2);
    g12 := Integers()!G[1][2]; g13 := Integers()!G[1][3]; g23 := Integers()!G[2][3];
    b1 := Integers()!Gg[1]; b2 := Integers()!Gg[2]; b3 := Integers()!Gg[3];
    c0 := Integers()!cst; Am := h33 mod pw;
    hist := [];
    for b in [0..pw-1] do
        row := [0 : v in [1..pw]];
        for t in [0..pw-1] do
            v := (Am*t*t + b*t) mod pw; row[v+1] +:= 1;
        end for;
        Append(~hist, row);
    end for;
    N := 0;
    for x1 in [0..pw-1] do for x2 in [0..pw-1] do
        B := (g13*x1 + g23*x2 - b3) mod pw;
        C := (h11*x1*x1 + h22*x2*x2 + g12*x1*x2 - b1*x1 - b2*x2 + c0) mod pw;
        N +:= hist[B+1][((-C) mod pw)+1];
    end for; end for;
    return N, pw;
end function;

// beta = 0 density product g_X(n) and the S_c-dependent en/ed, for a lattice with Gram GX
gpart := function(GX, n)
    detX := Determinant(ChangeRing(GX, Integers()));
    Sc := Sort([p : p in Set(PrimeDivisors(detX)) join Set(PrimeDivisors(Numerator(n)))]);
    D0 := -(Numerator(n)*Denominator(n));
    K := QuadraticField(D0); dd := Discriminant(Integers(K));
    chi := KroneckerCharacter(dd);
    g := Rationals()!1;
    z := Vector(Rationals(), [0,0,0]);
    for p in Sc do
        w := 1 + 2*Valuation(Rationals()!(2*n), p);
        cnt, pw := rep_count(GX, z, n, p, w);
        g *:= (Rationals()!cnt)/(Rationals()!pw)^2;
        if g eq 0 then break; end if;
    end for;
    en := &*[Rationals() | 1 - Evaluate(chi, p)/p : p in Sc];
    ed := &*[Rationals() | 1 - 1/(Rationals()!p)^2 : p in Sc];
    return g, en/ed, Sc;
end function;

emit("START");
for base in [<15,2>, <21,2>] do
    D, N := Explode(base);
    Ld := ShimuraCurveLattice(D, N);
    Qc := ChangeRing(Ld`Q, Integers()); Qr := ChangeRing(Ld`Q, Rationals());
    denom := Ld`denom; disc_grp := Ld`disc_grp; to_disc := Ld`to_disc;
    detL := Determinant(Qc);
    emit(Sprintf("=== X0^%o(%o), det L = %o ===", D, N, detL));

    // isotropic beta != 0 (all have N_beta = 2 on these bases) and their overlattice Grams
    overs := [];
    for eta in disc_grp do
        v := ChangeRing(eta@@to_disc, Rationals());
        if not IsIntegral((v*Qr, v)/(2*denom^2)) then continue; end if;
        w := v/denom;
        if IsZero(w) then continue; end if;
        Nb := LCM([Denominator(x) : x in Eltseq(w)]);
        if Nb ne 2 then emit(Sprintf("  SKIP beta with N_beta = %o (needs primitive-chi Thm 1.4)", Nb)); continue; end if;
        // L_beta = Z^3 + Z*beta, with 2*beta integral: 2*L_beta = HNF(2I, 2beta)
        rows := [Vector(Integers(), [2*x : x in Eltseq(Vector(Rationals(), [i eq j select 1 else 0 : j in [1..3]]))]) : i in [1..3]];
        M4 := Matrix(Integers(), 4, 3, [Eltseq(r) : r in rows] cat [[Integers()!(2*x) : x in Eltseq(w)]]);
        H := HermiteForm(M4);
        Hb := Matrix(Rationals(), 3, 3, [Eltseq(H[i]) : i in [1..3]])/2;   // basis of L_beta
        GB := Hb * Qr * Transpose(Hb);
        // must be an even integral lattice of index 2 (det scales by 1/4)
        ok := &and[IsIntegral(GB[i][j]) : i,j in [1..3]] and &and[IsIntegral(GB[i][i]/2) : i in [1..3]];
        emit(Sprintf("  beta = %-20o  det L_beta = %o (expect %o)  even-integral? %o",
                     Eltseq(w), Determinant(GB), detL/4, ok));
        if ok then Append(~overs, ChangeRing(GB, Integers())); end if;
    end for;
    emit(Sprintf("  #overlattices built = %o", #overs));

    // compare, per n, the beta=0-only coefficient vs the full isotropic-sum G coefficient
    emit("  n :  g_L (en/ed)_L  |  sum over overlattices  |  ratio  b^G_0(n)/b^L_0(n)");
    for n in [1..12] do
        gL, rL, ScL := gpart(Qc, Rationals()!n);
        bL := rL*gL;                                  // times the common 96/(D N) * cond * h/w
        bG := -bL;                                    // the  - E_{A,0}  term
        for GB in overs do
            gB, rB, ScB := gpart(GB, Rationals()!n);
            bG +:= 2*rB*gB;                           // factor 2 = sqrt(|det L|/|det L_beta|)
        end for;
        emit(Sprintf("  n=%-3o  b^L = %-12o  b^G = %-12o  %o", n, bL, bG,
                     bL eq 0 select (bG eq 0 select "both 0" else "L=0 but G!=0")
                              else Sprintf("ratio = %o", bG/bL)));
    end for;
end for;
emit("DONE");
