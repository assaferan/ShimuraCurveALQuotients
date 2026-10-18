// X_0^39(11)/W (genus one) by the degree-two recipe of Guo-Yang, Compositio 153 (2017),
// Example 42, as run for X_0^142(1)* in tests/GuoYangExample42.m.
//
// ⚠ NOT YET RUN.  The form ring here has level M = 4 * 429 = 1716 (four times the odd part of
// D*N), whose polytope solutions are not in polymake/ and whose first solve timed out on the Mac;
// the ring is being built on lovelace.  Run this from the repository root once the
// polymake_solution_1716_* files exist, with NORMALIZ_BIN exported:
//
//     magma -b tools/genus1_39_11.m < /dev/null > genus1_39_11.log 2>&1
//     magma -b stage:=0 tools/genus1_39_11.m < /dev/null        # object checks only, no form ring
//     magma -b stage:=1 tools/genus1_39_11.m < /dev/null        # through the feasibility sweep
//     magma -b nonfund:=1 tools/genus1_39_11.m < /dev/null      # also tabulate the point -28
//     magma -b D:=142 N:=1 tools/genus1_39_11.m < /dev/null    # dry run on Example 42 (a minute)
//
// The dry run on X_0^142(1)* exercises every stage against the known answer (E142A1; the forms,
// the table and the curve are checked independently in tests/GuoYangExample42.m) and was run on
// 2026-10-18; it is how the stages below were debugged without the level-1716 form ring.
//
// What is known about the object (checked 2026-10-18, stage 0 below re-checks it):
//   * the star quotient has genus 1 and X_0^39(11)/{1,11,39,429} has genus 3;
//   * the Jacobian is 429b1: y^2 + xy = x^3 - 24x + 63, rank 1, torsion of order 4;
//   * rational CM points (class number <= 8): -7, -19, -24, -28, -123, -187, -195, -403, -715,
//     -3003.  ⚠ -28 is NOT fundamental (conductor 2 over -7).  On main the point pool
//     (RationalandQuadraticCMPoints) lists fundamental discriminants only, because Schofer's
//     formula at a non-maximal order is not validated here (the non-fundamental corrections live
//     on the kappa0-proof branch), so the value table below uses the nine fundamental points and
//     admits -28 only with nonfund:=1, as a check rather than as input;
//   * the branch points of the genus-3 cover are the two conjugate pairs of discriminants -52 and
//     -132, each an elliptic point of order 2 on the star (fixed by w_13 and w_33); -52 has class
//     number 2 and -132 class number 4, so the search's point pool is taken up to class number 4.
//
// The torsion matters: a divisor may be principal only after multiplication by 2 or 4, in which
// case the rational search returns a form with fractional principal part whose Borcherds product
// is a root of a function, not a function.  IntegralSolution keeps only integral representatives.
//
// Stages.
//   0  object checks and the CM pool;
//   1  feasibility of every divisor of the shape Z(d1) + Z(d2) - 2 Z(d_oo) for d_oo rational and
//      d1, d2 rational, and of Z(-52) - 2 Z(d_oo), Z(-132) - 2 Z(d_oo), all in ONE call so the form
//      ring is built once (AllowMissing reports which divisors have no form);
//   2  from two poles that each carry two forms with disjoint zero sets, the coordinates x, y at
//      the rational CM points, the bidegree-(2, 2) relation and the Weierstrass model, which
//      must be 429b1;
//   3  the values at those points of the branch-pair forms Z(-52) - 2 Z(d_oo) and
//      Z(-132) - 2 Z(d_oo), ready to be expressed in x and y.
// Every stage prints what it found; nothing here asserts a value that has not been checked by
// hand, because none has been computed yet.

AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 2);
SetColumns(0);
if not assigned stage then stage := 3; else stage := StringToInteger(stage); end if;

if assigned D then D := StringToInteger(D); N := StringToInteger(N); else D := 39; N := 11; end if;
// the object: Cremona label of the Jacobian, the rational CM points to tabulate, the quadratic
// cycles whose divisors Z(d) - 2 Z(d_oo) are wanted (the branch points at 39_11)
if <D, N> eq <39, 11> then
    label := "429b1";
    table_ds := [-7, -19, -24, -123, -187, -195, -403, -715, -3003];
    if assigned nonfund then Insert(~table_ds, 4, -28); end if;     // see the header
    branch := [-52, -132];
    pool_bd := 4;
elif <D, N> eq <142, 1> then
    label := "142a1";
    table_ds := [-3, -4, -8, -19, -20, -24, -40, -43, -148, -232];
    branch := [-91];      // a quadratic cycle, to exercise the path; Z(-91) - 2 Z(d_oo) is not expected to be principal
    pool_bd := 2;
else
    error "no object data for this (D, N)";
end if;
W := Set(Divisors(D*N));
Xstar := CreateShimuraQuot(D, N, W);
Xstar`g := GenusShimuraCurveQuotient(D, N, W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};

// ---- stage 0: the object -------------------------------------------------------------------
printf "X_0^%o(%o)*: genus %o; genera of the quotients by the index-2 subgroups: %o\n", D, N, star`g,
       [<curves[i]`W, curves[i]`g> : i in star`CoveredBy];
assert star`g eq 1;
Ejac := EllipticCurve(CremonaDatabase(), label);
printf "%o: %o, rank %o, torsion %o\n", label, Ejac, Rank(Ejac), #TorsionSubgroup(Ejac);
rat, quad := RationalandQuadraticCMPoints(star : bd := pool_bd, coprime_to_level := false);
printf "rational CM points (class number <= %o): %o\nquadratic: %o\n", pool_bd, [p[1] : p in rat], [p[1] : p in quad];
ells := NumberOfEllipticPointsByCMOrder(star);
for q in Keys(ells) do printf "elliptic points of order %o: %o\n", q, [<d, ells[q][d]> : d in Keys(ells[q])]; end for;
for d in table_ds do
    assert DegreeOfFieldOfDefinitionOfCMPoint(star, d) eq 1;
end for;
rat_pool := [p[1] : p in rat];
assert Set(branch) subset {p[1] : p in quad};
if stage eq 0 then quit; end if;

// ---- stage 1: feasibility -------------------------------------------------------------------
// Only small discriminants are affordable as divisor support: the cusp-oo window must reach
// |m| = |d| (d odd) or |d|/4, and the cusp-0 ladder climbs by D0 * |m|.
small := [d : d in rat_pool | Abs(d) le 200];
poles := small;
pairs := [[small[i], small[j]] : j in [i+1..#small], i in [1..#small]];
divisors := [];
for P in poles do
    for pr in pairs do
        if P in pr then continue; end if;
        Append(~divisors, <pr, P>);
    end for;
    for b in branch do Append(~divisors, <[b], P>); end for;
end for;
printf "stage 1: %o candidate divisors over poles %o\n", #divisors, poles;
t0 := Realtime();
fs := BorcherdsForms(star, curves : Divisors := divisors, PoolClassNumber := pool_bd,
                                    IntegralSolution := true, AllowMissing := true);
printf "stage 1 done in %os: %o of %o divisors have a form\n", Realtime() - t0, #Keys(fs), #divisors;
found := Sort([k : k in Keys(fs)]);
for k in found do
    printf "  key %o: divisor %o   principal part at oo %o\n", k, divisors[k], qExpansionAtoo(fs[k], 1);
end for;
printf "  without a form: %o\n", [divisors[k] : k in [1..#divisors] | k notin Keys(fs)];
if stage eq 1 then quit; end if;

// ---- stage 2: coordinates, relation, curve -------------------------------------------------
// For each pole, the forms with two rational zeros; two with disjoint zero sets give a coordinate.
coord_choices := [];
for P in poles do
    ks := [k : k in found | divisors[k][2] eq P and #divisors[k][1] eq 2];
    disjoint := [<a, b> : a in ks, b in ks | a lt b and IsEmpty(Set(divisors[a][1]) meet Set(divisors[b][1]))];
    if not IsEmpty(disjoint) then Append(~coord_choices, <P, disjoint[1][1], disjoint[1][2]>); end if;
end for;
printf "stage 2: coordinates available at poles %o\n", [c[1] : c in coord_choices];
if #coord_choices lt 2 then
    printf "stage 2: fewer than two poles carry a coordinate; stopping.\n";
    quit;
end if;
// two coordinates with distinct poles generate the function field (a function in C(x) with a
// double pole elsewhere would have a pole at x's pole too)
cx := coord_choices[1]; cy := coord_choices[2];
t0 := Realtime();
vals := BorcherdsFormValues([fs[cx[2]], fs[cx[3]], fs[cy[2]], fs[cy[3]]], table_ds, star);
printf "stage 2: Schofer values in %os\n", Realtime() - t0;
for i in [1..4] do printf "  |psi_%o| at %o: %o\n", i, table_ds, vals[i]; end for;
xs, cx0, cx1 := CoordinateFromBorcherdsForms(vals[1], vals[2], table_ds, divisors[cx[2]][1], divisors[cx[3]][1]);
ys, cy0, cy1 := CoordinateFromBorcherdsForms(vals[3], vals[4], table_ds, divisors[cy[2]][1], divisors[cy[3]][1]);
printf "x (pole %o, zeros %o, ones %o; c_0 = %o, c_1 = %o):\n  %o\n", cx[1], divisors[cx[2]][1], divisors[cx[3]][1], cx0, cx1, xs;
printf "y (pole %o, zeros %o, ones %o; c_0 = %o, c_1 = %o):\n  %o\n", cy[1], divisors[cy[2]][1], divisors[cy[3]][1], cy0, cy1, ys;
F := BidegreeTwoRelation(xs, ys);
printf "relation: %o = 0\n", F;
E, toE := EllipticCurveFromRelation(F, xs, ys);
printf "minimal model: %o\n  Cremona %o (expected %o); isomorphic: %o\n", E, CremonaReference(E), label, IsIsomorphic(E, Ejac);
if stage eq 2 then quit; end if;

// ---- stage 3: the branch pairs ---------------------------------------------------------------
for k in found do
    if #divisors[k][1] eq 1 then
        v := BorcherdsFormValues([fs[k]], table_ds, star)[1];
        printf "branch-pair form %o: |psi| at %o: %o\n", divisors[k], table_ds, v;
    end if;
end for;
quit;
