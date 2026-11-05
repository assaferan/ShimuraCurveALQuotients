// An exact model of a star quotient X_0^D(N)/W of GENUS ONE, with its rational CM points placed
// on the curve, by the degree-two recipe of Guo-Yang, Compositio 153 (2017), Example 42
// (reproduced against the paper in tests/GuoYangExample42.m).
//
//     magma -b D:=158 N:=1 tools/genus1star.m < /dev/null > genus1star_158_1.log 2>&1
//
// Options (name:=value): stage (0 object checks; 1 the divisor search; 2 the model, default;
// 3 also the quadratic cycles that are fibres of a coordinate), bd (class-number bound of the CM
// pool, default 4), maxpole / maxzero (largest |d| used as a pole / a zero of a Borcherds form,
// defaults 100 / 200), intsol (1, default: keep only integral principal parts, so each Borcherds
// form is a function and not a root of one; 0 to accept rational ones), exclude (discriminants
// left out of the value table, comma-separated, e.g. exclude:=-51), verb (verbosity, default 1).
//
// What is printed, and how it is checked.
//   * The Jacobian is identified independently of the model: by Jacquet-Langlands the star's
//     Jacobian is the elliptic curve of the rational newform of level D*M (M | N) whose
//     Atkin-Lehner signs make it quaternionically invariant, i.e. a_p = +1 at p | D (the sign is
//     flipped at the discriminant primes) and a_p = -1 at p || N.  The model's curve must have
//     those traces of Frobenius at every good prime up to 100, which identifies the isogeny class;
//     the model itself identifies the curve in it, since the coordinates are functions on the star
//     (see the next item).
//   * Guo-Yang's Assumption 40 (header of GenusOneStar.m): it is not needed when E(Q) is
//     torsion-free, since the divisor class of a form with integral principal part has finite
//     order and is then trivial; it is a genuine hypothesis when E(Q) has rational torsion, and
//     the driver says which case it is in.  With rational 2-torsion and an Atkin-Lehner quotient
//     of genus one available, the driver can instead build each coordinate as psi^2, psi a form
//     with divisor A - P: that divisor is not principal, so its class is the unique nonzero point
//     of E(Q)[2] and psi^2 is a function.  That quotient is then provably the cover on which psi
//     becomes a function, so the signs it supplies are forced (SquareCoordinateSigns, and the
//     argument in the header of GenusOneStar.m); the model is unconditional in that case too.
//   * Every CM point of fundamental discriminant and class number <= bd that is rational on the
//     star (the pool of RationalandQuadraticCMPoints, elliptic points included) is placed on the
//     minimal model, as a point and as an element of the Mordell-Weil group (origin: the CM point
//     named as base point).  A point whose Schofer value is not a rational number stops the run
//     (BorcherdsFormValues); leaving it out with exclude:= makes the model conditional on the
//     values at the remaining points.  Schofer's formula at a discriminant not coprime to N has
//     no theoretical warrant here (see the note at CandidateDiscriminants).
//   * Stage 3 lists the quadratic CM cycles Z(d) = {P, P'} with Z(d) - 2 P_oo principal for some
//     rational CM point P_oo, i.e. P + P' = 2 P_oo exactly in the group: such a cycle is a fibre
//     of the coordinate with pole P_oo.  (The branch points of a double cover X/<w_m> -> X* are
//     the images of the fixed points of the other involutions, CM points by orders containing
//     sqrt(-m') with m' the other elements of W; tests of hyperellipticity need those, see
//     BiellipticBranchCriterion in GenusOneStar.m.)
//   * The Borcherds forms, the Schofer values and the relation are printed so that a value can be
//     traced; the file ends with one line  RESULT := <...>;  that Magma can read back.
//
// Limits.  Non-squarefree N is refused upstream (the CM-field code asserts).  On odd D the
// cusp-0 ladder reaches the class-number-2 pool plus the points used, so "no form" there may
// mean "not reached" rather than "not principal"; and each rung m of that ladder builds a pool
// of pole order D0*|m| at the zero cusp, whose size is what bounds the memory of a run
// (BorcherdsForms drops the rungs of unused pool points above MaxZeroPoleOrder, so keep maxzero
// small on odd D).

AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
if not assigned D or not assigned N then error "usage: magma -b D:=.. N:=.. tools/genus1star.m"; end if;
D := StringToInteger(D);
N := StringToInteger(N);
excluded := assigned exclude select [StringToInteger(x) : x in Split(exclude, ",")] else [];
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
// every rational CM point supports divisors and enters the value table unless excluded
support_ds := Sort([p[1] : p in rat], func<a, b | Abs(a) - Abs(b)>);
table_ds := [d : d in support_ds | d notin excluded];
if #excluded gt 0 then printf "left out of the value table on request: %o\n", excluded; end if;
// the characters of the Borcherds forms (header of GenusOneStar.m): a nontrivial one needs an
// etale double cover of the star among the index-2 quotients
// Assumption 40 is automatic when the Jacobian has no rational torsion.  Otherwise the square-
// coordinate route removes it, provided the torsion is Z/2 on every curve of the isogeny class --
// the star may be any of them -- and the star has a genus-one Atkin-Lehner quotient, which is then
// provably the cover on which a form of divisor A - P becomes a function (header of
// GenusOneStar.m).  Full rational 2-torsion would leave that identification open, so it is
// excluded here rather than only its exponent.
etale := EtaleDoubleCovers(star, curves);
// over every curve of every candidate isogeny class, since the star may be any of them
cls := IsEmpty(cands) select [] else &cat[IsogenousCurves(EllipticCurve(f)) : f in cands];
etors := IsEmpty(cls) select 1 else LCM([Exponent(TorsionSubgroup(C)) : C in cls]);
two_torsion := [];
for C in cls do
    T := TorsionSubgroup(C);
    Append(~two_torsion, #[Pt : Pt in T | 2*Pt eq T!0 and Pt ne T!0]);
end for;
printf "over the %o curve(s) of the candidate class(es): torsion exponent %o, nonzero rational 2-torsion points %o\n",
       #cls, etors, two_torsion;
single_two := not IsEmpty(two_torsion) and forall{t : t in two_torsion | t eq 1};
if etors eq 1 then
    printf "the Jacobian is torsion-free: every Borcherds form with integral principal part is a function here, and the model does not depend on Assumption 40\n";
elif etors eq 2 and single_two and not IsEmpty(etale) then
    printf "torsion Z/2 on every curve of the class, and the genus-one quotient(s) %o: the coordinates will be squares of forms with divisor A - P, whose signs that quotient then determines\n",
           [Sort(Setseq(c`W)) : c in etale];
    error if #etale gt 1, "several genus-one Atkin-Lehner quotients: the signs cannot be attributed to one of them";
else
    printf "⚠ the square route is not available here (torsion exponent %o, nonzero rational 2-torsion points %o, genus-one quotients %o): the model below is CONDITIONAL on Assumption 40\n",
           etors, two_torsion, #etale;
end if;
square_route := (etors eq 2) and single_two and not IsEmpty(etale);
sq := [];
if stage eq 0 then quit; end if;

// ---- stage 1: the divisor search ------------------------------------------------------------
poles := [d : d in support_ds | Abs(d) le maxpole];
small := [d : d in support_ds | Abs(d) le maxzero];
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
if square_route then
    // a form with divisor A - P exists exactly when A - P has finite order in E(Q); its square,
    // the function with divisor 2A - 2P, spans L(2P) with the constants
    deg1 := [<[A], P> : A in support_ds, P in poles | A ne P];
    t0 := Realtime();
    f1 := BorcherdsForms(star, curves : Divisors := deg1, PoolClassNumber := bd, IntegralSolution := true, AllowMissing := true);
    partners := [deg1[k] : k in Sort(Setseq(Keys(f1)))];
    printf "stage 1 (squares): degree-one divisors with a form, %os: %o\n", Realtime() - t0, partners;
    error if IsEmpty(partners), "no degree-one divisor has a form: no square-coordinate is available";
    divisors := [<[pr[1][1], pr[1][1]], pr[2]> : pr in partners];
    for P in poles do
        for q in quad do
            if Abs(q[1]) le maxzero then Append(~divisors, <[q[1]], P>); end if;
        end for;
    end for;
end if;
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
if square_route then
    // x = the function with divisor 2A - 2P, normalised to 1 at a table point, its sign at the
    // other points from the fields of definition on the etale cover (one global sign is free and
    // only changes x into -x)
    // every divisor the search is given has rational CM points as its zeros and its pole, which
    // the identification of the cover depends on
    sq := [k : k in found | #divisors[k][1] eq 2 and divisors[k][1][1] eq divisors[k][1][2] and divisors[k][2] in table_ds];
    printf "stage 2 (squares): square-coordinates at the poles %o\n", [divisors[k][2] : k in sq];
    error if #sq lt 2, "fewer than two poles carry a square-coordinate";
    cx := <divisors[sq[1]][2], sq[1], sq[1]>; cy := <divisors[sq[2]][2], sq[2], sq[2]>;
    t0 := Realtime();
    vals := BorcherdsFormValues([fs[sq[1]], fs[sq[2]]], table_ds, star);
    printf "stage 2: Schofer values in %os\n", Realtime() - t0;
    for i in [1..2] do printf "  |x_%o| at %o: %o\n", i, table_ds, vals[i]; end for;
    signs := SquareCoordinateSigns(etale[1], table_ds);
    printf "  signs from the fields of definition on the cover: %o\n", signs;
    coord := function(v)
        j0 := [j : j in [1..#table_ds] | Type(v[j]) ne Infty and v[j] ne 0][1];
        return [* Type(v[j]) eq Infty select Infinity() else signs[j]*signs[j0]*v[j]/v[j0] : j in [1..#table_ds] *];
    end function;
    xs := coord(vals[1]); ys := coord(vals[2]);
    printf "x (pole %o, double zero %o):\n  %o\ny (pole %o, double zero %o):\n  %o\n",
           cx[1], divisors[sq[1]][1][1], xs, cy[1], divisors[sq[2]][1][1], ys;
else
coord_choices := [];
for P in poles do
    // a coordinate needs its zeros in the value table (they fix its normalisation)
    ks := [k : k in found | divisors[k][2] eq P and #divisors[k][1] eq 2 and divisors[k][1] subset table_ds];
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
end if;
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
for f in matches do
    nc := #IsogenousCurves(EllipticCurve(f));
    if nc eq 1 then
        printf "  the class %o has one curve: isomorphic to it: %o\n", CremonaReference(EllipticCurve(f)), IsIsomorphic(E, EllipticCurve(f));
    else
        // the coordinates are functions on the star when the square route was taken or the
        // Jacobian is torsion-free, and then the model is the star curve itself; otherwise the
        // forms may be roots of functions and the curve is picked out only under Assumption 40
        printf "  the class of %o has %o curves; the model picks out %o%o\n",
               CremonaReference(EllipticCurve(f)), nc, CremonaReference(E),
               (square_route or etors eq 1) select "" else " (under Assumption 40)";
    end if;
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
        // a square-coordinate has 2P = 2A with P != A, so A itself is not the pole
        if z[1] eq z[2] then halves := [Q : Q in halves | Q ne pts[z[1]]]; end if;
        // distinct CM points of the star are distinct points of the model, so a candidate
        // already carrying another discriminant is not this pole
        taken := {pts[d] : d in Keys(pts) | d ne P};
        if #halves gt 1 then halves := [Q : Q in halves | Q notin taken]; end if;
        if #halves eq 1 then
            pts[P] := halves[1]; pos[P] := mw_coords(halves[1]); progress := true;
        else
            printf "  the pole %o is one of %o (2-torsion ambiguity); left unplaced\n", P, halves;
        end if;
    end for;
end while;
if square_route then
    // the partner A of each pole P must differ from P by a point of order 2 on the model
    for k in sq do
        A := divisors[k][1][1]; P := divisors[k][2];
        if IsDefined(pts, A) and IsDefined(pts, P) then
            T := pts[A] - pts[P];
            printf "  A - P for (A, P) = (%o, %o): %o, order %o\n", A, P, T, Order(T);
            assert Order(T) eq 2;
        end if;
    end for;
end if;
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
