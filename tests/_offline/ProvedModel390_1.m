// THE STORED MODEL OF X_0^390(1)/{1,6,26,39} INVOLVES NO CHOICE.  About 15 minutes.
//
// The curve is the fibre product, over the star X_0^390(1)*, of its two double covers above
// {1,5,6,26,30,39,130,195} and {1,2,3,6,13,26,39,78}, stored in data/models/models_390_1.m as
// y^2 = f(s) in one hauptmodul s.  Over the algebraic closure the curve depends only on the zeros of
// the two f, so the model is proved once each f is determined, up to a constant, by exact data.  The
// route (AllEquationsAboveCovers) has three inputs, and this test checks that none of them is a
// selection:
//   * the absolute values of s, 1 - s and y^2 at the CM points, from Schofer's formula, a theorem
//     at N = 1 where it has no m = 0 term;
//   * their signs: at every rational CM point exactly one sign pattern of s and 1 - s is
//     consistent, and exactly one sign of y^2 is a square in the field of definition of the point
//     on the cover; no sign is left to the search over sign combinations;
//   * the linear algebra: the constraints y^2(P) = f(s(P)) have a one-dimensional solution space,
//     so f is unique up to the constant (the twist, which hyperellipticity does not see).
// It also checks that the stored f are the ones these data give.

printf "the model of X_0^390(1)/{1,6,26,39} involves no choice...";
SetColumns(0);
P<x> := PolynomialRing(Rationals());
curves := GetHyperellipticCandidates();
Xstar := rep{X : X in curves | X`D eq 390 and X`N eq 1 and IsStarCurve(X)};
fs := GetBorcherdsForms(Xstar, curves : Prec := 100);
d_divs := &cat[[T[1] : T in DivisorOfBorcherdsForm(f, Xstar)] : f in [fs[-1], fs[-2]]];
all_cm_pts := CandidateDiscriminants(Xstar, curves : Keep := Set(d_divs));
num_vals := Maximum([2*curves[i]`g + 5 : i in Xstar`CoveredBy]);
abs_tab, all_cm_pts := AbsoluteValuesAtCMPoints(Xstar, curves, all_cm_pts, fs :
                         MaxNum := num_vals, Prec := 100, Exclude := {}, Include := Set(d_divs));
ReduceTable(abs_tab);
tab := ValuesAtCMPoints(abs_tab, all_cm_pts);
assert #tab`AmbiguousSigns eq 0 and IsEmpty(tab`UnscaledKeys);
ds := tab`Discs; V := tab`Values; si := tab`sIndex; ti := tab`sTildeIndex;
flds := abs_tab`FldsOfDefn;
// the hauptmodul: s + s~ = 1, and at each rational point one sign pattern only
for j in [1..#ds] do
    v := V[si][j]; w := V[ti][j];
    if Type(v) eq Infty then assert Type(w) eq Infty; continue; end if;
    assert Type(v) eq FldRatElt and v + w eq 1;
    if v eq 0 or w eq 0 then continue; end if;      // s = 0 or 1: a zero of a Borcherds product, no sign
    assert #[<e1, e2> : e1, e2 in [-1, 1] | e1*Abs(v) + e2*Abs(w) eq 1] eq 1;
end for;
kernels := RationalConstraintsOnEquations(tab, curves);
crvs, _, keys, _ := EquationsOfCovers(tab, all_cm_pts);
models := eval (Read("data/models/models_390_1.m") cat "\nreturn models;");
for W in [{Integers()|1,5,6,26,30,39,130,195}, {Integers()|1,2,3,6,13,26,39,78}] do
    i := [r : r in [1..#tab`Keys_fs] | tab`Keys_fs[r] gt 0 and curves[tab`Keys_fs[r]]`W eq W][1];
    for j in [1..#ds] do
        v := V[i][j];
        if Type(v) ne FldRatElt or v eq 0 then continue; end if;
        F := flds[tab`Keys_fs[i]][ds[j]];
        assert #[e : e in [-1, 1] | IsSquare(F!(e*v))] eq 1;
    end for;
    assert #kernels[Index(tab`K_idxs, i)] eq 1;
    f := HyperellipticPolynomials(crvs[[m : m->k in keys | curves[k]`W eq W][1]]);
    stored := models[Sort(Setseq(W))][1][2];
    c := LeadingCoefficient(stored) / LeadingCoefficient(f);
    assert stored eq c*f and IsSquare(c);           // the same cover, up to the twist by a square
end for;
printf " ok\n";
