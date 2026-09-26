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
// PART 3 checks every candidate against the oracle over a range of levels, D > 1 included.
// PART 4 runs the ModSym check end to end on two [FH] Theorem 4 curves, pinned to the involution
// recorded for them in data/curves_after_UpdateCurves8.dat (the first decisive candidate).

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
assert nControl eq 5;

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
    for W in ALSubgroups(L) do
        V_names, idx_sets, names, good := ModularNonALInvolutionCandidates(D, N, W);
        mats := [(v eq "S2") select S2 else ((v eq "V2") select get_V2(L) else get_V3(L)) : v in V_names];
        for i->I in idx_sets do
            v := &*[mats[j] : j in I];
            for o in good[i] do
                assert isinv(v*al_matrix(o, L), W, L);
                nChecked +:= 1;
            end for;
        end for;
    end for;
end for;
assert nChecked gt 500;

// ------------------------------------------------ PART 4: end to end, [FH] Theorem 4
for t in [<63, {1, 9}, "V3 W1">, <120, {1, 5, 24, 120}, "V2 W1">] do
    X := CreateShimuraQuot(1, t[1], t[2]);
    X`g := GenusShimuraCurveQuotient(1, t[1], t[2]);
    r, nm := CheckModularNonALInvolutionModSym(X);
    assert r eq 1 and nm eq t[3];
end for;

printf "done (%o candidates checked).\n", nChecked;
