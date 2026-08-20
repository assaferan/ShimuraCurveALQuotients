// TEST THE INCOHERENT HYPOTHESIS, properly.
//
// V represents a over Q_p  <=>  g = V perp <-a> (rank 4) is ISOTROPIC over Q_p
//                          <=>  det(g) is not a square, OR det(g) is a square and eps_p(g) = (-1,-1)_p.
// This rank-4 criterion avoids the ternary Hasse-invariant conventions.  (NOTE: there is no parity
// constraint on |Diff_finite| here -- the even-cardinality theorem needs det(V perp <-a>) to be a
// GLOBAL square, which it is not in general.  An earlier version of this script asserted parity and
// was itself wrong.)
//
// TWO READINGS OF THE INCOHERENT PREDICTION, both tested:
//   C_oo positive definite (3,0): it represents m > 0, so oo is not in Diff and |Diff| = 1 means
//        EXACTLY ONE finite p fails;
//   C_oo negative definite (0,3): it fails on m > 0, so oo is in Diff and |Diff| = 1 means
//        NO finite p fails.
AttachSpec("ShimuraQuotients.spec");

function is_square_Qp(a, p)
    v := Valuation(a, p);
    if IsOdd(v) then return false; end if;
    u := a / p^v;
    num := Numerator(u); den := Denominator(u);
    if p eq 2 then return ((num * InverseMod(den, 8)) mod 8) eq 1; end if;
    return IsSquare(GF(p) ! (num * InverseMod(den, p)));
end function;

function represents_Qp(diag3, a, p)
    diag := diag3 cat [-a];
    d := &*diag;
    if not is_square_Qp(d, p) then return true; end if;
    eps := &*[Integers() | HilbertSymbol(Rationals()!diag[i], Rationals()!diag[j], p)
              : i in [1..4], j in [1..4] | i lt j];
    return eps eq HilbertSymbol(Rationals()!(-1), Rationals()!(-1), p);
end function;

btrue := AssociativeArray();
btrue[[15,2]] := [<2, 4>, <10, -4>, <30, 0>];
btrue[[6,5]]  := [<10, -6>, <15, -2>, <30, -6>];
btrue[[10,3]] := [<3, -2>, <12, -2>, <30, -6>];

for b in [[15,2],[6,5],[10,3]] do
    D := b[1]; N := b[2];
    Ld := ShimuraCurveLattice(D, N);
    G := ChangeRing(Ld`Q, Rationals());
    Dg := OrthogonalizeGram(G);
    diag3 := [Dg[i][i] : i in [1..3]];
    printf "\n===== X0^%o(%o)  diag(Q) = %o, det = %o =====\n", D, N, diag3, &*diag3;
    // AUDIT: at a prime not dividing 2*det the ternary form is unimodular, hence isotropic, hence
    // represents everything.  This is a valid check on the criterion (unlike parity).
    allp := Sort(Setseq({p : p in PrimeDivisors(2*3*5*7*11*13*Integers()!(&*diag3))}));
    for p in allp do
        if (2*Integers()!(&*diag3)) mod p eq 0 then continue; end if;
        for a in [1..12] do
            error if not represents_Qp(diag3, Rationals()!a, p),
                Sprintf("good prime %o wrongly fails to represent %o", p, a);
        end for;
    end for;
    printf "  good-prime audit passed\n";
    printf "  %-5o %-7o %-16o %-16o %o\n", "r", "b_true", "Diff_fin(Q)", "Diff_fin(-Q)", "A: |Diff|=1  B: |Diff|=0";
    for t in btrue[b] do
        r := Rationals()!t[1];
        dQ  := [p : p in allp | not represents_Qp(diag3, r, p)];
        dmQ := [p : p in allp | not represents_Qp([-x : x in diag3], r, p)];
        // incoherent prediction: b nonzero iff Diff_finite is EMPTY (so Diff = {oo} alone)
        okA := ((#dQ eq 1) and (t[2] ne 0)) or ((#dQ ne 1) and (t[2] eq 0));
        okB := ((#dQ eq 0) and (t[2] ne 0)) or ((#dQ ne 0) and (t[2] eq 0));
        printf "  %-5o %-7o %-16o %-16o A:%-6o B:%-6o\n", t[1], t[2], dQ, dmQ,
               okA select "ok" else "FAIL", okB select "ok" else "FAIL";
    end for;
end for;
quit;
