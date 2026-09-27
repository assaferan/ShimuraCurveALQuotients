// Degrees of the fields of definition of CM points on the Atkin-Lehner quotient
// X^D_0(N)/W, read off Pic(R) and the Atkin-Lehner combinatorics of W alone (no
// ring class field), and the rational and degree-2 CM points built from them.
// For the fields themselves see FieldsOfDefinitionOfCMPointFast in SchoferFormula.m.
//
// The degree of the field of definition of the CM points by R on X^D_0(N)/W is
//
//     2 * h(R) * 2^eps / #W_gal / (2 if complex conjugation is active else 1)
//
// where W_gal <= W is the subgroup acting as Galois automorphisms (Lemma 5.9 of
// Gonzalez-Rotger), 2^eps is the order of the part of W_gal acting trivially --
// the Atkin-Lehner involutions that fix R itself -- and the last factor records
// whether some element of W identifies a point with its complex conjugate.
// For the star quotient of a squarefree N with D = 1 the conjugation factor is
// always present and this reduces to the familiar h(R)/2^(omega(N/N_R) - eps).


intrinsic CMClassLists() -> Assoc
{The imaginary quadratic discriminants of each class number up to 8, keyed by
 class number: every order, maximal or not, and every class group, not only the
 2-groups (a point on a proper subquotient can be rational without an elementary
 abelian 2 class group).}
    CNs := AssociativeArray();
    CNs[1] := {-3,-4,-7,-8,-11,-12,-16,-19,-27,-28,-43,-67,-163};
    CNs[2] := {-15,-20,-24,-32,-35,-36,-40,-48,-51,-52,-60,-64,-72,-75,-88,-91,-99,-100,-112,-115,-123,-147,-148,-187,-232,-235,-267,-403,-427};
    CNs[3] := {-23,-31,-44,-59,-76,-83,-92,-107,-108,-124,-139,-172,-211,-243,-268,-283,-307,-331,-379,-499,-547,-643,-652,-883,-907};
    CNs[4] := {-39,-55,-56,-63,-68,-80,-84,-96,-120,-128,-132,-136,-144,-155,-156,-160,-168,-171,-180,-184,-192,-195,-196,-203,-208,-219,-220,-228,-240,-252,-256,-259,-275,-280,-288,-291,-292,-312,-315,-323,-328,-340,-352,-355,-363,-372,-387,-388,-400,-408,-435,-448,-475,-483,-507,-520,-532,-555,-568,-592,-595,-603,-627,-667,-708,-715,-723,-760,-763,-772,-795,-928,-955,-1003,-1012,-1027,-1227,-1243,-1387,-1411,-1435,-1467,-1507,-1555};
    CNs[5] := {-47,-79,-103,-127,-131,-179,-188,-227,-316,-347,-412,-443,-508,-523,-571,-619,-683,-691,-739,-787,-947,-1051,-1123,-1723,-1747,-1867,-2203,-2347,-2683};
    CNs[6] := {-87,-104,-116,-135,-140,-152,-175,-176,-200,-204,-207,-212,-216,-244,-247,-300,-304,-324,-339,-348,-364,-368,-396,-411,-424,-432,-436,-451,-459,-460,-472,-484,-492,-496,-515,-531,-540,-588,-628,-648,-675,-676,-688,-700,-707,-747,-748,-771,-808,-828,-835,-843,-856,-867,-891,-931,-940,-963,-988,-1048,-1059,-1068,-1072,-1075,-1083,-1099,-1107,-1108,-1147,-1192,-1203,-1219,-1267,-1315,-1323,-1347,-1363,-1432,-1563,-1588,-1603,-1612,-1675,-1708,-1843,-1915,-1963,-2227,-2283,-2403,-2443,-2515,-2563,-2608,-2787,-2923,-3235,-3427,-3523,-3763,-4075};
    CNs[7] := {-71,-151,-223,-251,-284,-343,-463,-467,-487,-587,-604,-811,-827,-859,-892,-1163,-1171,-1372,-1483,-1523,-1627,-1787,-1852,-1948,-1987,-2011,-2083,-2179,-2251,-2467,-2707,-3019,-3067,-3187,-3907,-4603,-5107,-5923};
    CNs[8] := {-95,-111,-164,-183,-224,-248,-260,-264,-272,-276,-295,-299,-308,-320,-336,-360,-371,-376,-380,-384,-392,-395,-420,-444,-452,-456,-468,-480,-504,-512,-528,-539,-544,-548,-552,-564,-576,-579,-580,-583,-600,-612,-616,-624,-632,-640,-651,-660,-672,-712,-720,-732,-736,-768,-784,-792,-819,-820,-832,-840,-852,-868,-880,-900,-904,-912,-915,-939,-952,-960,-979,-987,-995,-1008,-1024,-1032,-1035,-1043,-1060,-1092,-1120,-1128,-1131,-1152,-1155,-1156,-1168,-1180,-1195,-1204,-1240,-1248,-1252,-1275,-1288,-1299,-1312,-1320,-1332,-1339,-1348,-1360,-1380,-1395,-1408,-1428,-1443,-1488,-1528,-1540,-1552,-1587,-1600,-1632,-1635,-1651,-1659,-1672,-1683,-1731,-1752,-1768,-1771,-1780,-1792,-1795,-1803,-1827,-1828,-1848,-1864,-1912,-1939,-1947,-1992,-1995,-2020,-2035,-2059,-2067,-2080,-2088,-2107,-2115,-2128,-2139,-2163,-2212,-2248,-2272,-2275,-2307,-2308,-2323,-2332,-2368,-2392,-2395,-2419,-2451,-2475,-2587,-2611,-2632,-2667,-2715,-2755,-2788,-2827,-2832,-2907,-2947,-2968,-2995,-3003,-3040,-3088,-3172,-3243,-3283,-3315,-3355,-3403,-3448,-3507,-3595,-3627,-3712,-3787,-3843,-3883,-3963,-4048,-4123,-4195,-4267,-4323,-4387,-4747,-4843,-4867,-5083,-5467,-5587,-5707,-5947,-6307,-7987};
    return CNs;
end intrinsic;

// Discriminants of the orders fixed by the Atkin-Lehner involution w_m: the
// fixed points of w_m are CM by Z[sqrt(-m)] (discriminant -4m) and, when
// -m = 1 mod 4, by Z[(1+sqrt(-m))/2] (discriminant -m).  This depends only on
// m, not on how m divides D versus N, so it applies unchanged on a Shimura
// curve.  w_2 additionally fixes the points with CM by Z[i].
function ALFixedDiscs(m)
    if m eq 1 then return {}; end if;
    if m eq 2 then return {-4, -8}; end if;
    if (m mod 4) eq 3 then return {-m, -4*m}; end if;
    return {-4*m};
end function;

intrinsic DegreeOfFieldOfDefinitionOfCMPoint(X::ShimuraQuot, d::RngIntElt) -> RngIntElt
{Degree over Q of the field of definition of the CM points by the order of
 discriminant d on X, or 0 if X carries no CM point by that order.  Computed
 from Pic(R) and the Atkin-Lehner data alone, without building a ring class
 field; agrees with Degree of the fields returned by FieldsOfDefinitionOfCMPoint.}
    D := X`D;
    N := X`N;
    W := X`W;
    R := QuadraticOrder(BinaryQuadraticForms(d));
    f := Conductor(R);
    chi := KroneckerCharacter(d);

    D_R      := &*[Integers()| p : p in PrimeDivisors(D) | chi(p) eq -1];
    N_R      := &*[Integers()| p : p in PrimeDivisors(N) | chi(p) eq 1 or (f mod p eq 0)];
    N_star_R := &*[Integers()| p : p in PrimeDivisors(N) | chi(p) eq 1 and (f mod p ne 0)];

    // Proposition 5.6 plus the GCD(D, f) = 1 correction: an order non-maximal at
    // a prime ramified in the quaternion algebra has embedding number 0, which
    // the congruence alone does not see.
    if ((Discriminant(R) mod ((D*N) div (D_R*N_star_R))) ne 0) or (GCD(D, f) ne 1) then
        return 0;
    end if;

    // W_gal: the elements of W that act as Galois automorphisms (Lemma 5.9).
    W_gal := {mm : mm in W | ((D*N) div (D_R*N_R)) mod mm eq 0};

    // The involutions of W_gal that fix R act trivially on the ring class field,
    // so they do not cut the field of definition down.  They form a subgroup;
    // its order is the 2-power correction 2^eps.
    fixers := {mm : mm in W_gal | d in ALFixedDiscs(mm)} join {1};
    ker := AllALsFromGens(fixers, D*N) meet W_gal;
    eps_factor := #ker;

    // Complex conjugation is realised by w_m composed with an Artin symbol, for
    // m = D_R * N_star_R.  It halves the degree exactly when some mm in W has
    // m * mm Galois; m * mm need not itself lie in W (e.g. W = {1, D*N}).
    m := D_R*N_star_R;
    cc_active := exists{mm : mm in W |
        ((D*N) div (D_R*N_R)) mod AtkinLehnerMul(m, mm, D*N) eq 0};

    h_R := #PicardGroup(R);
    // Exact: #W_gal / eps_factor is the order of a subgroup of Gal(H_R/K).
    // Assert it so that integer division cannot hide a wrong W_gal or eps.
    assert (2 * h_R * eps_factor) mod #W_gal eq 0;
    deg := 2 * h_R * eps_factor div #W_gal;
    if cc_active then
        assert deg mod 2 eq 0;
        deg div:= 2;
    end if;
    return deg;
end intrinsic;

// Number of CM points by R lying over one point of the quotient: the 2^omega
// Atkin-Lehner orbit, shortened on the involutions that fix R.
function CMOrbitDenominator(X, d)
    W_size := #X`W;
    fixers := {mm : mm in X`W | d in ALFixedDiscs(mm)} join {1};
    return W_size div #(AllALsFromGens(fixers, X`D*X`N) meet X`W);
end function;

intrinsic RationalCMDiscs(X::ShimuraQuot : max_class_num := 0) -> Assoc
{Discriminants of the orders whose CM points on X are rational, mapped to the
 number of such points.  Uses Shimura reciprocity through
 DegreeOfFieldOfDefinitionOfCMPoint; requires D*N squarefree.  max_class_num > 0
 considers only class numbers up to that bound.}
    require IsSquarefree(X`D * X`N) : "RationalCMDiscs requires D*N squarefree";
    CNs := CMClassLists();
    cm_pts := AssociativeArray();
    hs := Sort(Setseq(Keys(CNs)));
    if max_class_num gt 0 then
        hs := [h : h in hs | h le max_class_num];
    end if;
    for h in hs do
        for disc in CNs[h] do
            R := QuadraticOrder(BinaryQuadraticForms(disc));
            opt_emb := NumberOfOptimalEmbeddings(R, X`D, X`N);
            if opt_emb eq 0 then continue; end if;
            if DegreeOfFieldOfDefinitionOfCMPoint(X, disc) ne 1 then continue; end if;
            vprintf ShimuraQuotients, 2: "\tfound rational CM point at discriminant %o\n", disc;
            cm_pts[disc] := opt_emb / CMOrbitDenominator(X, disc);
        end for;
    end for;
    return cm_pts;
end intrinsic;

intrinsic Degree2Points(X::ShimuraQuot : max_class_num := 0, compute_fields := true) -> Assoc
{Discriminants of the orders whose CM points on X have degree 2, mapped to
 <multiplicity, fields, field discriminants>.  max_class_num > 0 considers only
 class numbers up to that bound.  compute_fields := false leaves the two field
 slots as false, skipping the ring class field build that dominates the cost.}
    require IsSquarefree(X`D * X`N) : "Degree2Points requires D*N squarefree";
    CNs := CMClassLists();
    cm_pts := AssociativeArray();
    hs := Sort(Setseq(Keys(CNs)));
    if max_class_num gt 0 then
        hs := [h : h in hs | h le max_class_num];
    end if;
    for h in hs do
        for disc in CNs[h] do
            R := QuadraticOrder(BinaryQuadraticForms(disc));
            opt_emb := NumberOfOptimalEmbeddings(R, X`D, X`N);
            if opt_emb eq 0 then continue; end if;
            if DegreeOfFieldOfDefinitionOfCMPoint(X, disc) ne 2 then continue; end if;
            mult := opt_emb / CMOrbitDenominator(X, disc);
            if not compute_fields then
                cm_pts[disc] := <mult, false, false>;
                continue;
            end if;
            flds := FieldsOfDefinitionOfCMPointFast(X, disc);
            error if #flds eq 0,
                Sprintf("degree 2 at discriminant %o but no field of definition (D=%o, N=%o)",
                        disc, X`D, X`N);
            cm_pts[disc] := <mult, flds, [* Discriminant(F) : F in flds *]>;
        end for;
    end for;
    return cm_pts;
end intrinsic;

intrinsic CountDegree2Points(X::ShimuraQuot) -> FldRatElt
{Total number of degree-2 CM points on X.}
    cm_pts := Degree2Points(X : compute_fields := false);
    return &+[Rationals()| cm_pts[k][1] : k in Keys(cm_pts)];
end intrinsic;
