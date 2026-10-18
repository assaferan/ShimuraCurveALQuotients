// An exact model of a star quotient X_0^D(N)/W of GENUS ONE, with its rational CM points placed
// on the curve, by the degree-two recipe of Guo-Yang, Compositio 153 (2017), Example 42
// (reproduced against the paper in tests/GuoYangExample42.m).
//
//     magma -b D:=158 N:=1 tools/genus1star.m < /dev/null > genus1star_158_1.log 2>&1
//
// Options (name:=value): stage (0 object checks; 1 the divisor search; 2 the model, default;
// 3 also the quadratic branch pairs), bd (class-number bound of the CM pool, default 4),
// maxpole / maxzero (largest |d| used as a pole / a zero of a Borcherds form, defaults 100 / 200),
// intsol (1, default: keep only integral principal parts, so each Borcherds form is a function and
// not a root of one; 0 to accept rational ones), verb (verbosity, default 1).
//
// What is printed, and how it is checked.
//   * The Jacobian is identified independently of the model: by Jacquet-Langlands the star's
//     Jacobian is the elliptic curve of the rational newform of level D*M (M | N) whose
//     Atkin-Lehner signs make it quaternionically invariant, i.e. a_p = +1 at p | D (the sign is
//     flipped at the discriminant primes) and a_p = -1 at p || N.  The model's curve must have
//     those traces of Frobenius at every good prime up to 100, and when the isogeny class has a
//     single curve it must be isomorphic to it.
//   * Every rational CM point of class number <= bd is placed on the minimal model, as a point
//     and as an element of the Mordell-Weil group (origin: the CM point named as base point).
//   * Stage 3 lists the quadratic CM cycles Z(d) = {P, P'} with Z(d) - 2 P_oo principal for some
//     rational CM point P_oo, which says P + P' = 2 P_oo exactly in the group.  These are the
//     branch points of the genus-three double covers of the star (tests of hyperellipticity use
//     them: see BiellipticBranchCriterion in GenusOneStar.m).
//   * The Borcherds forms, the Schofer values and the relation are printed so that a value can be
//     traced; the file ends with one line  RESULT := <...>;  that Magma can read back.
//
// Known limits.  Non-squarefree N is refused upstream (the CM-field code asserts); on odd D the
// cusp-0 ladder only reaches the class-number-2 pool plus the points used, so "no form" there may
// mean "not reached" rather than "not principal" (tools history: 15_2).  Schofer's formula at a
// discriminant not coprime to N has no theoretical warrant here, only the record of the N > 1
// bases that passed their external checks (see the note at CandidateDiscriminants).
//
// X_0^39(11)* (Jacobian 429b1, branch discriminants -52 and -132 of the genus-3 quotient by
// {1,11,39,429}) and X_0^55(7)* are the targets this was written for; their form rings (levels
// 1716 and 1540) are being built on lovelace and are not to be attempted on the Mac.

AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
if not assigned D then error "usage: magma -b D:=.. N:=.. tools/genus1star.m"; end if;
D := StringToInteger(D); N := StringToInteger(N);
stage := assigned stage select StringToInteger(stage) else 2;
bd := assigned bd select StringToInteger(bd) else 4;
maxpole := assigned maxpole select StringToInteger(maxpole) else 100;
maxzero := assigned maxzero select StringToInteger(maxzero) else 200;
intsol := assigned intsol select StringToInteger(intsol) eq 1 else true;
SetVerbose("ShimuraQuotients", assigned verb select StringToInteger(verb) else 1);

W := {m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1};
Xstar := CreateShimuraQuot(D, N, W);
Xstar`g := GenusShimuraCurveQuotient(D, N, W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};

// ---- stage 0: the object -------------------------------------------------------------------
printf "X_0^%o(%o)*: genus %o; quotients by the index-2 subgroups: %o\n", D, N, star`g,
       [<Sort(Setseq(curves[i]`W)), curves[i]`g> : i in star`CoveredBy];
error if star`g ne 1, "the star quotient does not have genus one";
cands := [];
for M in Divisors(N) do
    for orb in Newforms(CuspForms(Gamma0(D*M), 2)) do
        f := orb[1];
        if Degree(BaseRing(Parent(f))) ne 1 then continue; end if;
        if not forall{p : p in PrimeDivisors(D) | Coefficient(f, p) eq 1} then continue; end if;
        if not forall{p : p in PrimeDivisors(M) | Valuation(D*M, p) gt 1 or Coefficient(f, p) eq -1} then continue; end if;
        Append(~cands, f);
    end for;
end for;
printf "Jacobian candidates by Atkin-Lehner signs: %o\n",
       [<Level(f), CremonaReference(EllipticCurve(f)), #IsogenousCurves(EllipticCurve(f))> : f in cands];
rat, quad := RationalandQuadraticCMPoints(star : bd := bd, coprime_to_level := false);
ells := NumberOfEllipticPointsByCMOrder(star);
printf "rational CM points (class number <= %o): %o\nquadratic: %o\nelliptic points: %o\n",
       bd, [p[1] : p in rat], [p[1] : p in quad],
       [<q, [<d, ells[q][d]> : d in Keys(ells[q])]> : q in Keys(ells) | #Keys(ells[q]) gt 0];
table_ds := Sort([p[1] : p in rat], func<a, b | Abs(a) - Abs(b)>);
if stage eq 0 then quit; end if;

// ---- stage 1: the divisor search ------------------------------------------------------------
poles := [d : d in table_ds | Abs(d) le maxpole];
small := [d : d in table_ds | Abs(d) le maxzero];
pairs := [[small[i], small[j]] : j in [i+1..#small], i in [1..#small]];
divisors := [];
for P in poles do
    for pr in pairs do
        if P notin pr then Append(~divisors, <pr, P>); end if;
    end for;
    for q in quad do
        if Abs(q[1]) le maxzero then Append(~divisors, <[q[1]], P>); end if;
    end for;
end for;
printf "stage 1: %o candidate divisors over the poles %o\n", #divisors, poles;
t0 := Realtime();
fs := BorcherdsForms(star, curves : Divisors := divisors, PoolClassNumber := bd,
                                    IntegralSolution := intsol, AllowMissing := true);
printf "stage 1 done in %os: %o of %o divisors have a form\n", Realtime() - t0, #Keys(fs), #divisors;
found := Sort([k : k in Keys(fs)]);
for k in found do
    printf "  key %o: divisor %o   principal part at oo %o\n", k, divisors[k], qExpansionAtoo(fs[k], 1);
end for;
if stage eq 1 then quit; end if;

// ---- stage 2: coordinates, relation, curve, CM points -----------------------------------------
coord_choices := [];
for P in poles do
    ks := [k : k in found | divisors[k][2] eq P and #divisors[k][1] eq 2];
    disjoint := [<a, b> : a in ks, b in ks | a lt b and IsEmpty(Set(divisors[a][1]) meet Set(divisors[b][1]))];
    if not IsEmpty(disjoint) then Append(~coord_choices, <P, disjoint[1][1], disjoint[1][2]>); end if;
end for;
printf "stage 2: coordinates available at the poles %o\n", [c[1] : c in coord_choices];
error if #coord_choices lt 2, "fewer than two poles carry a coordinate";
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
// base point: the first table point that is finite on the plane model and nonsingular there
E, toE := EllipticCurveFromRelation(F, xs, ys);
Cp := Domain(toE);
finite := [j : j in [1..#table_ds] | Type(xs[j]) ne Infty and Type(ys[j]) ne Infty];
base_j := [j : j in finite | IsNonsingular(Cp, Cp![xs[j], ys[j], 1])][1];
printf "minimal model: %o\n  Cremona %o; base point (origin) the CM point of discriminant %o\n", E, CremonaReference(E), table_ds[base_j];

// the Jacobian check, independent of the model: traces of Frobenius against the newform
good_ps := [p : p in PrimesUpTo(100) | (2*D*N) mod p ne 0];
matches := [f : f in cands | forall{p : p in good_ps | TraceOfFrobenius(E, p) eq Coefficient(f, p)}];
printf "Jacobian check: traces at %o good primes match the newform(s) %o\n", #good_ps,
       [CremonaReference(EllipticCurve(f)) : f in matches];
error if IsEmpty(matches), "the model's curve matches no Jacobian candidate";
single := [f : f in matches | #IsogenousCurves(EllipticCurve(f)) eq 1];
for f in single do
    printf "  the class %o has one curve: isomorphic to it: %o\n", CremonaReference(EllipticCurve(f)), IsIsomorphic(E, EllipticCurve(f));
end for;

// the CM points on E: the finite ones through the map, a pole of x or of y through its divisor
// (x has divisor zeros - 2 pole, so the pole is half the sum of the two zeros, up to 2-torsion)
G, mw := MordellWeilGroup(E);
printf "E(Q) = %o\n", G;
// coordinates in G by a bounded search (the map has no preimage rule); CM points of small
// discriminant have small height, so a box of 40 in each free direction is ample here
gens := [mw(G.i) : i in [1..Ngens(G)]];
ords := [Order(G.i) : i in [1..Ngens(G)]];
box := CartesianProduct([ords[i] eq 0 select [-40..40] else [0..ords[i]-1] : i in [1..Ngens(G)]]);
mw_coords := function(Pt)
    for n in box do
        if &+[E | n[i]*gens[i] : i in [1..#gens]] eq Pt then return [n[i] : i in [1..#gens]]; end if;
    end for;
    return "outside the search box";
end function;
pos := AssociativeArray();
pts := AssociativeArray();
for j in finite do
    pts[table_ds[j]] := toE(Cp![xs[j], ys[j], 1]);
    pos[table_ds[j]] := mw_coords(pts[table_ds[j]]);
end for;
// a pole P of a found form with two placed zeros z1, z2 satisfies 2P = z1 + z2; repeat while
// it places something, since the pole of x may be a zero of the form placing the pole of y
progress := true;
while progress do
    progress := false;
    for k in found do
        P := divisors[k][2]; z := divisors[k][1];
        if IsDefined(pts, P) or #z ne 2 or not (IsDefined(pts, z[1]) and IsDefined(pts, z[2])) then continue; end if;
        halves := [Q : Q in DivisionPoints(pts[z[1]] + pts[z[2]], 2)];
        if #halves eq 1 then
            pts[P] := halves[1]; pos[P] := mw_coords(halves[1]); progress := true;
        else
            printf "  the pole %o is one of %o (2-torsion ambiguity); left unplaced\n", P, halves;
        end if;
    end for;
end while;
printf "CM points on E (discriminant, point, Mordell-Weil coordinates):\n";
for d in table_ds do
    if IsDefined(pts, d) then printf "  %o: %o   %o\n", d, pts[d], pos[d]; end if;
end for;
placed := [<d, Eltseq(pts[d]), pos[d]> : d in table_ds | IsDefined(pts, d)];
if stage eq 2 then
    printf "RESULT := <%o, %o, %m, \"%o\", %o, %m>;\n", D, N, aInvariants(E), CremonaReference(E), table_ds[base_j], placed;
    quit;
end if;

// ---- stage 3: quadratic cycles that are x-fibres --------------------------------------------
branch := [];
for k in found do
    if #divisors[k][1] eq 1 then
        q := divisors[k][1][1]; P := divisors[k][2];
        printf "Z(%o) - 2 P_%o is principal: the two conjugate points of discriminant %o sum to 2 P_%o%o\n",
               q, P, q, P, IsDefined(pts, P) select Sprintf(" = %o", 2*pts[P]) else "";
        Append(~branch, <q, P>);
    end if;
end for;
printf "RESULT := <%o, %o, %m, \"%o\", %o, %m, %m>;\n", D, N, aInvariants(E), CremonaReference(E), table_ds[base_j], placed, branch;
quit;
