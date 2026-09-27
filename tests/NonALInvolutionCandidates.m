// tests/NonALInvolutionCandidates.m -- which v*W_o the non-AL filter is allowed to test.
//
// CheckModularNonALInvolution{ModSym,Trace} read the genus g' of Y/<v W_o>, Y = X_0(D,N)/W, as
// the quotient by an INVOLUTION: g' = 0 => hyperelliptic, and fix = 2g - 4g' + 2 fixed points.
// Both readings are false for an element of order 3 or 4, so the candidate list must contain
// only involutions of Y.  Before this test the list included S2 V2 (never an involution with
// the o it was paired with) and V2 V3 W_o for the wrong o, and it skipped genuine involutions
// such as V3 W10 on X_0(90) (a per-prime test used where [FH] Lemma 1's eps is multiplicative).
//
// PART 1 checks the oracle IsModularInvolutionOnQuotient on hand-derived cases, including
// non-involutions it must reject (so a vacuous oracle cannot pass PART 2).  PART 2 pins the
// candidate lists at small levels to the values derived by hand from [FH] Lemma 1:
//   eps(d) = 1 iff the 3-free part of d is 2 mod 3;  (V3 W_o)^2 = W_9^eps(o);
//   (V2 V3 W_o)^2 = W_9^eps(2^a o), 2^a || N;  S2 W_o is an involution iff o is odd.
// PART 2b pins why the triple S2 V2 V3 is not listed although it gives involutions: they are
// conjugate by w_{2^a} to listed S2 V3 W_o'.
// PART 3 checks every candidate against the oracle over a range of levels, D > 1 included.
// PART 4 runs the ModSym check end to end on two [FH] Theorem 4 curves.  It pins the verdict and
// the quotient genus of the reported witness (recomputed by the trace formula), not the witness
// name: which decisive candidate is reported first is not part of the mathematics.

printf "NonALInvolutionCandidates.m: non-AL candidates are exactly the involutions...";

M2Z := MatrixAlgebra(Integers(), 2);
S2 := M2Z![2,1,0,2];
isinv := IsModularInvolutionOnQuotient;

// ------------------------------------------------ PART 1: the oracle, both directions
nControl := 0;
// order 3: S2 W_4 on X_0(4)
assert not isinv(S2*al_matrix(4, 4), {1}, 4);                        nControl +:= 1;
// order 4: S2 W_8 on X_0(8)
assert not isinv(S2*al_matrix(8, 8), {1}, 8);                        nControl +:= 1;
// S2 V2 on X_0(8): order 4
assert not isinv(S2*get_V2(8), {1}, 8);                              nControl +:= 1;
// V2 V3 on X_0(72) with 9 notin W: (V2 V3)^2 = W_9^eps(8) = W_9
assert not isinv(get_V2(72)*get_V3(72), {1}, 72);                    nControl +:= 1;
// V3 W_2 on X_0(18): (V3 W_2)^2 = W_9^eps(2) = W_9
assert not isinv(get_V3(18)*al_matrix(2, 18), {1}, 18);              nControl +:= 1;
// ... and each becomes an involution once W_9 is divided out, or with the right o
assert isinv(get_V2(72)*get_V3(72), {1, 9}, 72);
assert isinv(get_V3(18)*al_matrix(2, 18), {1, 9}, 18);
assert isinv(get_V3(90)*al_matrix(10, 90), {1}, 90);    // eps(10) = eps(2) eps(5) = 0
assert isinv(get_V2(72)*get_V3(72)*al_matrix(8, 72), {1}, 72);       // eps(8*8) = 0
assert isinv(get_V2(144)*get_V3(144), {1}, 144);        // 2^4 || 144: eps(16) = 0
assert isinv(S2*al_matrix(9, 36), {1}, 36);
// W_4 S2 W_4^-1 = [1,0;2,1]: an involution of X_0(4) that fails only the lower-left mod-L
// test of Gamma_0(4) membership, so this pins that part of the oracle.
assert isinv(M2Z![1,0,2,1], {1}, 4);
// Non-automorphisms that square into Gamma_0(L) and fix W = {1} trivially.  The first two were
// accepted before the determinant and Gamma_0(L)-conjugation checks; each is caught by both.
assert not isinv(M2Z![1,0,1,-1], {1}, 36);          // det -1
assert not isinv(M2Z![-6,-5,-3,-3], {1}, 9);        // det 3, 3 is not Q * square for Q in {1, 9}
// ... and one caught by each check alone: diag(1,-1) normalizes Gamma_0(L) but has det -1 (it
// does not preserve the upper half-plane); S = [0,-1;1,0] has det 1 but does not normalize Gamma_0(9).
assert not isinv(M2Z![1,0,0,-1], {1}, 36);
assert not isinv(M2Z![0,-1,1,0], {1}, 9);
nControl +:= 4;
assert nControl eq 9;

// ------------------------------------------------ PART 2: pinned candidate lists
V_names, idx_sets, names, good := ModularNonALInvolutionCandidates(1, 72, {1});
assert V_names eq ["S2", "V2", "V3"];
assert names eq ["S2", "V2", "V3", "S2 V3", "V2 V3"];                 // no S2 V2, each pair once
assert good eq [{1, 9}, {1, 8, 9, 72}, {1, 9}, {1, 9}, {8, 72}];

_, _, names, good := ModularNonALInvolutionCandidates(1, 72, {1, 9});
assert names eq ["S2", "V2", "V3", "S2 V3", "V2 V3"];
assert good eq [{1}, {1, 8, 72}, {1, 8, 72}, {1}, {1, 8, 72}];     // 9 in W: no eps condition

_, _, names, good := ModularNonALInvolutionCandidates(1, 90, {1});
assert names eq ["V3"];
assert good eq [{1, 9, 10, 90}];                                      // 10, 90 were skipped before

_, _, names, good := ModularNonALInvolutionCandidates(10, 9, {1});    // D > 1: level 90
assert names eq ["V3"];
assert good eq [{1, 9, 10, 90}];

// V3 does not descend to X_0(18)/<w_2>, and S2 does not descend when W has an even element
_, _, names, _ := ModularNonALInvolutionCandidates(1, 18, {1, 2});
assert names eq [];
_, _, names, _ := ModularNonALInvolutionCandidates(1, 36, {1, 4});
assert names eq ["V3"];

// ------------------------------------------------ PART 2b: S2 V2 V3 is omitted, and why
// S2 V2 V3 W_o IS an involution of X_0(72) for o in {8, 72}, but w_8 conjugates it to the listed
// S2 V3 W_{o/8}, so the quotients are isomorphic.  Pin: the involution, the conjugacy modulo
// Q^x Gamma_0(72) W, and equality of the two quotient genera by the trace formula.
function InQGammaW(M, W, L)
    M2Q := MatrixAlgebra(Rationals(), 2);
    for w in W do
        X := (M2Q!M) * (M2Q!al_matrix(w, L))^-1;
        ok, c := IsSquare(Determinant(X));
        if Determinant(X) gt 0 and ok and &and[IsIntegral(x) : x in Eltseq(X / c)]
           and (Integers()!(X / c)[2,1]) mod L eq 0 then return true; end if;
    end for;
    return false;
end function;
M2Q := MatrixAlgebra(Rationals(), 2);
T72 := M2Q!ModularInvolution("S2 V2 V3", 72);
SV72 := M2Q!ModularInvolution("S2 V3", 72);
W8 := M2Q!al_matrix(8, 72);
_, _, names, good := ModularNonALInvolutionCandidates(1, 72, {1});
assert "S2 V2 V3" notin names;
for o in [8, 72] do
    g := T72 * M2Q!al_matrix(o, 72);
    assert isinv(g, {1}, 72);                                         // a genuine involution ...
    h := SV72 * M2Q!al_matrix(o div 8, 72);
    assert (o div 8) in good[Index(names, "S2 V3")];                  // ... whose conjugate is listed
    assert InQGammaW(W8 * g * W8^-1 * h^-1, {1}, 72);                 // w_8 g w_8^-1 = S2 V3 W_{o/8}
    assert not InQGammaW(g * h^-1, {1}, 72);                          // (not equal without w_8)
    assert TraceDNewQuotient(Matrix(Integers(), T72), "S2 V2 V3", o, {1}, 1, 72)
        eq TraceDNewQuotient(Matrix(Integers(), SV72), "S2 V3", o div 8, {1}, 1, 72);
end for;
// and it is never an involution for o with 8 notdivides o
assert forall{o : o in [1, 9] | not isinv(T72 * M2Q!al_matrix(o, 72), {1}, 72)};

// ------------------------------------------------ PART 3: every candidate is an involution
function ALSubgroups(L)
    als := [Q : Q in Divisors(L) | GCD(Q, L div Q) eq 1];
    subs := {{1}};
    repeat
        n := #subs;
        for S in subs, a in als do
            Include(~subs, S join {AtkinLehnerMul(a, s, L) : s in S});
        end for;
    until #subs eq n;
    return subs;
end function;
nChecked := 0;
for DN in [<1,8>, <1,16>, <1,36>, <1,72>, <1,144>, <1,90>, <1,288>, <10,9>, <15,8>, <35,36>] do
    D := DN[1]; N := DN[2]; L := D*N;
    G0gens := Gamma0GeneratorMatrices(L);
    for W in ALSubgroups(L) do
        V_names, idx_sets, names, good := ModularNonALInvolutionCandidates(D, N, W);
        mats := [(v eq "S2") select S2 else ((v eq "V2") select get_V2(L) else get_V3(L)) : v in V_names];
        for i->I in idx_sets do
            v := &*[mats[j] : j in I];
            for o in good[i] do
                assert isinv(v*al_matrix(o, L), W, L : Gamma0Gens := G0gens);
                nChecked +:= 1;
            end for;
        end for;
    end for;
end for;
assert nChecked gt 500;

// ------------------------------------------------ PART 4: end to end, [FH] Theorem 4
// Genus of Y/<v W_o> for the witness named "v W_o", by the trace formula (independent of ModSym),
// after checking the named element is an involution of Y.
function WitnessQuotientGenus(nm, D, N, W)
    parts := Split(nm, " ");
    vname := &cat[(k eq 1 select "" else " ") cat parts[k] : k in [1..#parts-1]];
    o := StringToInteger(parts[#parts][2..#parts[#parts]]);
    V := ModularInvolution(vname, D*N);
    assert isinv(V*al_matrix(o, D*N), W, D*N);
    return TraceDNewQuotient(V, vname, o, W, D, N);
end function;
for t in [<63, {1, 9}>, <120, {1, 5, 24, 120}>] do
    X := CreateShimuraQuot(1, t[1], t[2]);
    X`g := GenusShimuraCurveQuotient(1, t[1], t[2]);
    r, nm := CheckModularNonALInvolutionModSym(X);
    assert r eq 1;
    gq := WitnessQuotientGenus(nm, 1, t[1], t[2]);
    assert (gq eq 0) or ((X`g eq 3) and (gq eq 2));
    r_tr := CheckModularNonALInvolutionTrace(X);
    assert r_tr eq 1;
end for;

printf "done (%o candidates checked).\n", nChecked;
