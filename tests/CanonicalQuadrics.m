// tests/CanonicalQuadrics.m
//
// The stage FilterByCanonicalQuadrics decides hyperellipticity by Max Noether's count: the pairwise
// products of a basis of holomorphic differentials of a curve of genus g >= 3 satisfy (g-1)(g-2)/2
// linear relations when the curve is hyperelliptic and (g-2)(g-3)/2 otherwise.
//
// [1] D = 1, against Ogg's list of the hyperelliptic modular curves X_0(N) (Ogg, "Hyperelliptic
//     modular curves", Bull. SMF 102 (1974), Theorem 2): among genus >= 3, X_0(N) is hyperelliptic
//     exactly for N = 30, 33, 35, 39, 40, 41, 48 (genus 3), 46, 47 (genus 4), 59 (genus 5), 71
//     (genus 6); X_0(34), (45), (64) (genus 3), (53) (genus 4) and (55) (genus 5) are not.  The
//     W-fixed count on the W-trivial quotient must reproduce this, and agree with Magma's own
//     IsHyperelliptic on the model from X0NQuotient.  N = 39 and 71 have index not divisible by 3.
// [2] D = 1 with W nontrivial: X_0(114)/<w_2, w_38>, X_0(130)/<w_2, w_26>, X_0(195)/<w_5, w_195>
//     are hyperelliptic of genus 3 by Magma's IsHyperelliptic on the model from X0NQuotient
//     (computed here); the stage must say so.
// [3] D > 1 from stored models: X_0(6,23)/w_23 is stored as a fibre product and is a smooth plane
//     quartic (tools/planequartic.m, PROVENANCE.md 2026-10-03); X_0(10,53)/w_530 is stored as
//     y^2 = f(x) of genus 3, accepted only after its genus and point counts are checked, and its
//     label says the model rests on the m = 0 term; a curve with no model file is left undecided;
//     a decided curve is not touched.
//
// The records are built directly (CreateShimuraQuot and the genus from the genus formula), with
// IsSubhyp unassigned, as the pipeline's undecided curves are.

printf "Testing the canonical-quadrics stage...";

// [1] Ogg's list, and Magma's IsHyperelliptic on X0NQuotient
for N in [30, 33, 35, 39, 40, 41, 46, 47, 48, 59, 71] do
    hyp, nrel, g := HyperellipticByQExpansions(N, {Integers()|1});
    assert hyp and g eq Genus(Gamma0(N)) and nrel eq (g-1)*(g-2) div 2;
    assert IsHyperelliptic(X0NQuotient(N, [Integers()|]));
end for;
for N in [34, 45, 64, 53, 55] do
    hyp, nrel, g := HyperellipticByQExpansions(N, {Integers()|1});
    assert (not hyp) and g eq Genus(Gamma0(N)) and nrel eq (g-2)*(g-3) div 2;
    assert not IsHyperelliptic(X0NQuotient(N, [Integers()|]));
end for;
// [2] the three quotients, by Magma's model
for c in [<114, [2, 38]>, <130, [2, 26]>, <195, [5, 195]>] do
    C := X0NQuotient(c[1], c[2]);
    assert Genus(C) eq 3 and IsHyperelliptic(C);
end for;

// records for the stage
function rec(D, N, W)
    X := CreateShimuraQuot(D, N, W);
    X`g := GenusShimuraCurveQuotient(D, N, W);
    return X;
end function;
curves := [rec(1, 40, {Integers()|1}), rec(1, 34, {Integers()|1}),
           rec(1, 114, {Integers()|1, 2, 19, 38}), rec(1, 130, {Integers()|1, 2, 13, 26}), rec(1, 195, {Integers()|1, 5, 39, 195}),
           rec(6, 23, {Integers()|1, 23}), rec(10, 53, {Integers()|1, 530}),
           rec(6, 109, {Integers()|1, 2}),            // no model file for this base
           rec(1, 48, {Integers()|1})];
assert [X`g : X in curves] eq [3, 3, 3, 3, 3, 3, 3, GenusShimuraCurveQuotient(6, 109, {Integers()|1, 2}), 3];
curves[9]`IsSubhyp := false; curves[9]`IsHyp := false; curves[9]`TestInWhichProved := "planted";   // a decided curve is left alone
FilterByCanonicalQuadrics(~curves);
expected := [true, false, true, true, true, false, true];
for i in [1..7] do
    assert assigned curves[i]`IsSubhyp and curves[i]`IsSubhyp eq expected[i] and curves[i]`IsHyp eq expected[i];
    assert "CanonicalQuadrics" in curves[i]`TestInWhichProved;
end for;
assert not assigned curves[8]`IsSubhyp;
assert curves[9]`TestInWhichProved eq "planted" and not curves[9]`IsSubhyp;
assert "stored hyperelliptic model" in curves[7]`TestInWhichProved and "m = 0 term" in curves[7]`TestInWhichProved;
assert "W-fixed cusp forms" in curves[1]`TestInWhichProved;
// a stored hyperelliptic model of the wrong genus is not accepted: the same curve against a file
// holding y^2 = x^5 + 1 (genus 2) under its key
dir := GetTempDir() cat "/cq_wrong_genus";
System("mkdir -p " cat dir);
P<x> := PolynomialRing(Rationals());
Write(dir cat "/models_10_53.m", "P<x> := PolynomialRing(Rationals());\nmodels := AssociativeArray();\nmodels[[1, 530]] := [ <2, P![1, 0, 0, 0, 0, 1], P![]> ];\n" : Overwrite := true);
wrong := [rec(10, 53, {Integers()|1, 530})];
FilterByCanonicalQuadrics(~wrong : ModelDir := dir);
assert not assigned wrong[1]`IsSubhyp;
printf " ok (Ogg's list and X0NQuotient, the three hyperelliptic quotients, two stored models, a wrong-genus model refused, one base without a file)\n";
