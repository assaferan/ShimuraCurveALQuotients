// Verify POINTWISE (per prime, per eta) whether the Schwagenscheidt Thm 1.4 local density
//    N_{gamma,n}(p^{w_p}) / p^{2 w_p},   w_p = 1 + 2 v_p(2 N_gamma n)
// equals LocalWhittakerAtOne(n, p, gamma, L, -Q). Earlier runs only compared the PRODUCT over S_c,
// so pointwise agreement must be confirmed before any test asserts it.

AttachSpec("ShimuraQuotients.spec");
fh := Open("pointwise_out.txt", "w");
emit := procedure(s) fprintf fh, "%o\n", s; Flush(fh); end procedure;

rep_count := function(G, gam, n, p, w)
    pw := p^w;
    if w le 0 then return 1, pw; end if;
    Gr := ChangeRing(G, Rationals());
    Gg := gam*Gr; cst := (gam*Gr, gam)/2 + n;
    if not IsIntegral(cst) then return -1, pw; end if;   // n not in Z - Q(gamma): undefined
    for i in [1..3] do if not IsIntegral(Gg[i]) then return -1, pw; end if; end for;
    h11 := Integers()!(G[1][1]/2); h22 := Integers()!(G[2][2]/2); h33 := Integers()!(G[3][3]/2);
    g12 := Integers()!G[1][2]; g13 := Integers()!G[1][3]; g23 := Integers()!G[2][3];
    b1 := Integers()!Gg[1]; b2 := Integers()!Gg[2]; b3 := Integers()!Gg[3];
    c0 := Integers()!cst; Am := h33 mod pw;
    hist := [];
    for b in [0..pw-1] do
        row := [0 : v in [1..pw]];
        for t in [0..pw-1] do v := (Am*t*t + b*t) mod pw; row[v+1] +:= 1; end for;
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

emit("START");
for base in [<15,2>, <21,2>] do
    D, N := Explode(base);
    Ld := ShimuraCurveLattice(D, N);
    Qint := ChangeRing(Ld`Q, Integers()); negQ := -Qint;
    Qr := ChangeRing(Ld`Q, Rationals()); denom := Ld`denom;
    Lfull := RSpaceWithBasis(IdentityMatrix(Integers(), 3));
    detprimes := Sort([p : p in Set(PrimeDivisors(Determinant(Qint)))]);
    emit(Sprintf("=== X0^%o(%o), det primes %o ===", D, N, detprimes));
    agree := 0; differ := 0; skipped := 0;
    z := Vector(Rationals(), [0,0,0]);
    // gamma = 0 (the oo-block case), integral n
    for n in [1..14] do
        for p in Sort([q : q in Set(detprimes) join Set(PrimeDivisors(n))]) do
            w := 1 + 2*Valuation(Rationals()!(2*n), p);
            cnt, pw := rep_count(Qint, z, Rationals()!n, p, w);
            if cnt lt 0 then skipped +:= 1; continue; end if;
            dens := (Rationals()!cnt)/(Rationals()!pw)^2;
            wh := LocalWhittakerAtOne(Rationals()!n, p, z, Lfull, negQ);
            if dens eq wh then
                agree +:= 1;
            else
                differ +:= 1;
                emit(Sprintf("  DIFFER n=%-3o p=%-3o w=%-2o : density=%-10o whittaker=%-10o ratio=%o",
                             n, p, w, dens, wh, wh eq 0 select "inf" else dens/wh));
            end if;
        end for;
    end for;
    emit(Sprintf("  gamma=0: agree=%o differ=%o skipped=%o", agree, differ, skipped));
end for;
emit("DONE");
