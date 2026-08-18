// Gap-1 probe: replace the local factor in the m=0 multiplier by the Schwagenscheidt
// (arXiv:1803.10550) Theorem 1.4 local factor, specialized to our case.
//
// Specialization (m = rank L = 3, k = 3/2, beta = 0, chi = 1 so N_beta = 1, g = 1):
//   eps = 1, c^finite = 1, and since X = p^{1-m/2-k} = p^-2 while p^{m-1} = p^2, the factor
//   (1 - p^{m-1} X) VANISHES, so
//        L^{(p)}_{gamma,n}(p^-2) = N_{gamma,n}(p^{w_p}) / p^{2 w_p},   w_p = 1 + 2 v_p(2 N_gamma n)
//   a pure local representation density (k = m/2 is the Siegel-Weil point). Here
//        N_{gamma,n}(a) = #{ r in L/aL : Q(r-gamma) + n = 0 mod a }.
// The Gram used is G = -Q_code (signature (2,1), det < 0), which is what makes the paper's
// discriminant D = 2(-1)^{(m+1)/2} N_gamma^2 n det(L) NEGATIVE, i.e. chi_{-n} -- matching the
// code's existing D0 = -(num*den) convention and the negQ already passed to LocalWhittakerAtOne.
//
// This probe prints, per (eta, n) term, the OLD local factor (LocalWhittakerAtOne) against the
// NEW density, and the resulting multiplier both ways.

AttachSpec("ShimuraQuotients.spec");
fh := Open("schwagen_out.txt", "w");
emit := procedure(s) fprintf fh, "%o\n", s; Flush(fh); end procedure;

// ---- N_{gamma,n}(p^w) via a histogram over the last coordinate -------------------------------
// F(x) = Q_G(x) - (x,gamma)_G + (Q_G(gamma) + n) is INTEGER-valued on Z^3 (L even, gamma in L').
// For fixed (x1,x2): F = A*x3^2 + B*x3 + C with A = G33/2 constant. We histogram
// t -> A t^2 + b t over t mod p^w for each b, then look up (-C).
rep_count := function(G, gam, n, p, w)
    pw := p^w;
    if w le 0 then return 1, pw; end if;
    Gg := gam*ChangeRing(G, Rationals());              // in Z^3 since gamma in L'
    cst := (gam*ChangeRing(G,Rationals()), gam)/2 + n; // Q_G(gamma) + n, in Z
    assert IsIntegral(cst);
    for i in [1..3] do assert IsIntegral(Gg[i]); end for;
    h11 := Integers()!(G[1][1]/2); h22 := Integers()!(G[2][2]/2); h33 := Integers()!(G[3][3]/2);
    g12 := Integers()!G[1][2]; g13 := Integers()!G[1][3]; g23 := Integers()!G[2][3];
    b1 := Integers()!Gg[1]; b2 := Integers()!Gg[2]; b3 := Integers()!Gg[3];
    c0 := Integers()!cst;
    Am := h33 mod pw;
    // hist[b+1][v+1] = #{ t mod pw : Am t^2 + b t = v mod pw }
    hist := [];
    for b in [0..pw-1] do
        row := [0 : v in [1..pw]];
        for t in [0..pw-1] do
            v := (Am*t*t + b*t) mod pw;
            row[v+1] +:= 1;
        end for;
        Append(~hist, row);
    end for;
    N := 0;
    for x1 in [0..pw-1] do
        for x2 in [0..pw-1] do
            B := (g13*x1 + g23*x2 - b3) mod pw;
            C := (h11*x1*x1 + h22*x2*x2 + g12*x1*x2 - b1*x1 - b2*x2 + c0) mod pw;
            N +:= hist[B+1][((-C) mod pw)+1];
        end for;
    end for;
    return N, pw;
end function;

// brute-force cross-check of rep_count (small moduli only)
rep_count_brute := function(G, gam, n, p, w)
    pw := p^w;
    Gr := ChangeRing(G, Rationals());
    Gg := gam*Gr; cst := (gam*Gr, gam)/2 + n;
    cnt := 0;
    for x1 in [0..pw-1] do for x2 in [0..pw-1] do for x3 in [0..pw-1] do
        x := Vector(Rationals(), [x1,x2,x3]);
        F := (x*Gr, x)/2 - (x*Gr, gam) + cst;
        if (Integers()!F) mod pw eq 0 then cnt +:= 1; end if;
    end for; end for; end for;
    return cnt;
end function;

// ---- the multiplier, with a switchable local factor ------------------------------------------
multiplier := function(foo, f0, Ldata, D, N, use_new)
    Q := Ldata`Q; disc_grp := Ldata`disc_grp; to_disc := Ldata`to_disc; denom := Ldata`denom;
    M := IsOdd(D*N) select 4*D*N else 2*D*N;
    Qint := ChangeRing(Q, Integers()); Qr := ChangeRing(Q, Rationals());
    negQ := -Qint;
    detprimes := Set(PrimeDivisors(Determinant(Qint)));
    Lfull := RSpaceWithBasis(IdentityMatrix(Integers(), 3));

    mod_M_to_vecs := AssociativeArray([0..M-1]);
    for j in [0..M-1] do mod_M_to_vecs[j] := []; end for;
    i0 := 0;
    for eta in disc_grp do
        if IsZero(eta) then i0 := eta; end if;
        v := ChangeRing(eta@@to_disc, Rationals());
        res := M*((v*Qr, v)/(2*denom^2));
        if not IsIntegral(res) then continue; end if;
        Append(~mod_M_to_vecs[Integers()!res mod M], eta);
    end for;

    contrib := function(eta, r, c)
        w_eta := ChangeRing(eta@@to_disc, Rationals())/denom;
        D0 := -(Numerator(r)*Denominator(r));
        K := QuadraticField(D0); dd := Discriminant(Integers(K));
        chi := KroneckerCharacter(dd);
        h := ClassNumber(K); wr := #TorsionSubgroup(UnitGroup(K));
        is_sq, cond_half := IsSquare(Rationals()!(r/AbsoluteValue(dd))); assert is_sq;
        Sc := Sort([p : p in detprimes join Set(PrimeDivisors(Numerator(r)))]);
        g := Rationals()!1;
        Ngam := (eta eq i0) select 1 else Order(eta);
        for p in Sc do
            if use_new ne 0 then
                // use_new = 1: Eisenstein lattice Gram -Q_code ; = 2: Gram +Q_code
                Guse := (use_new eq 1) select negQ else Qint;
                w := 1 + 2*Valuation(Rationals()!(2*Ngam*r), p);
                cnt, pw := rep_count(Guse, Vector(Rationals(), Eltseq(w_eta)), r, p, w);
                loc := (Rationals()!cnt)/(Rationals()!pw)^2;
            else
                loc := LocalWhittakerAtOne(r, p, Vector(Rationals(), Eltseq(w_eta)), Lfull, negQ);
            end if;
            emit(Sprintf("      p=%-3o  loc = %o", p, loc));
            g *:= loc;
        end for;
        en := &*[Rationals() | 1 - Evaluate(chi, p)/p : p in Sc];
        ed := &*[Rationals() | 1 - 1/(Rationals()!p)^2 : p in Sc];
        val := c * cond_half * (Rationals()!h/wr) * (en/ed) * g;
        emit(Sprintf("    term(N_gam=%o, n=%o, c=%o): h=%o w=%o cond=%o g=%o -> %o",
                     Ngam, r, c, h, wr, cond_half, g, val));
        return val;
    end function;

    T := Rationals()!0;
    for mm in [1..-Valuation(foo)] do
        c := Coefficient(foo, -mm);
        if c ne 0 then T +:= contrib(i0, Rationals()!mm, Rationals()!c); end if;
    end for;
    for j in [1..-Valuation(f0)] do
        c := Coefficient(f0, -j);
        if c eq 0 then continue; end if;
        r := (Rationals()!j)/M;
        for eta in mod_M_to_vecs[j mod M] do T +:= contrib(eta, r, Rationals()!c); end for;
    end for;
    return -96 * T / (D*N);
end function;

// ---- run ------------------------------------------------------------------------------------
emit("START");
R<q> := LaurentSeriesRing(Rationals());

// sanity: histogram count vs brute force on a small case
Ld15 := ShimuraCurveLattice(15, 2);
G15 := -ChangeRing(Ld15`Q, Integers());
z := Vector(Rationals(), [0,0,0]);
for pr in [<2,3>, <3,1>, <5,1>] do
    p, w := Explode(pr);
    a := rep_count(G15, z, Rationals()!2, p, w);
    b := rep_count_brute(G15, z, Rationals()!2, p, w);
    emit(Sprintf("selfcheck p=%o w=%o: hist=%o brute=%o %o", p, w, a, b, a eq b select "OK" else "MISMATCH"));
end for;

// X0^15(2) star hauptmodul principal part (as in tests/M0Multiplier.m); validated answer = 4
foo := 2*q^-10 - 2*q^-2 + 2*q^-1;
f0 := R!1;
emit("=== 15_2 OLD (LocalWhittakerAtOne) ===");
old15 := multiplier(foo, f0, Ld15, 15, 2, 0);
emit(Sprintf("15_2 OLD multiplier = %o   (validated target 4)", old15));
emit("=== 15_2 NEW, Eisenstein Gram = -Q_code ===");
new15a := multiplier(foo, f0, Ld15, 15, 2, 1);
emit(Sprintf("15_2 NEW(-Q) multiplier = %o", new15a));
emit("=== 15_2 NEW, Eisenstein Gram = +Q_code ===");
new15b := multiplier(foo, f0, Ld15, 15, 2, 2);
emit(Sprintf("15_2 NEW(+Q) multiplier = %o", new15b));
emit(Sprintf("SUMMARY: old=%o  new(-Q)=%o  new(+Q)=%o  target 4", old15, new15a, new15b));
emit("DONE");
