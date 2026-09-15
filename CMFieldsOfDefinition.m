// Fields of definition of the CM points on the Atkin-Lehner quotient
// X^D_0(N)/W, computed from ring class field data: their degree, and hence
// which discriminants give a point that is rational or of degree 2.
//
// Entry points:  DegreeOfFieldOfDefinitionOfCMPoint(X, d)   just the degree
//                RationalCMDiscs(X)                         discriminants of degree 1
//                Degree2Points(X)                           discriminants of degree 2
//                CMClassLists() / CMClassListsDeg2()        the candidate discriminants
//
// The degree is read off group orders alone, with no number field and no ring
// class field: only Pic(R) and the Atkin-Lehner combinatorics of W.  This is
// what makes scanning a whole class-number table affordable.  For the fields
// themselves see FieldsOfDefinitionOfCMPoint(Fast) in SchoferFormula.m, which
// computes the same degree the expensive way.
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
{Absolute values of the discriminants of the imaginary quadratic orders of each
 class number, keyed by class number.  Covers non-maximal orders as well as
 maximal ones, so it is the right candidate list for CM points.}
    CNs := AssociativeArray();
    CNs[1] := {3,4,7,8,11,12,16,19,27,28,43,67,163};
    CNs[2] := {15,20,24,32,35,36,40,48,51,52,60,64,72,75,88,91,99,100,112,115,123,147,148,187,232,
               235,267,403,427};
    CNs[4] := {84,96,120,132,160,168,180,192,195,228,240,280,288,312,315,340,352,372,408,435,448,
               483,520,532,555,595,627,708,715,760,795,928,1012,1435};
    CNs[8] := {420,480,660,672,840,960,1092,1120,1155,1248,1320,1380,1428,1540,1632,1848,1995,2080,
               3003,3040,3315};
    CNs[16] := {3360,5280,5460,7392};
    CNs[32] := {};
    CNs[64] := {};
    return CNs;
end intrinsic;

intrinsic CMClassListsDeg2() -> Assoc
{As CMClassLists, extended to the larger class numbers a degree-2 point can
 have.  Separate from CMClassLists because the degree-2 scan reaches one power
 of two further and the extra discriminants are not wanted in the rational scan.}
    CNs := AssociativeArray();
    CNs[1] := {3,4,7,8,11,12,16,19,27,28,43,67,163};
    CNs[2] := {15,20,24,32,35,36,40,48,51,52,60,64,72,75,88,91,99,100,112,115,123,147,148,187,232,
               235,267,403,427};
    CNs[4] := {39,55,56,63,68,80,84,96,120,128,132,136,144,155,156,160,168,171,180,184,192,195,196,
               203,208,219,220,228,240,252,256,259,275,280,288,291,292,312,315,323,328,340,352,355,
               363,372,387,388,400,408,435,448,475,483,507,520,532,555,568,592,595,603,627,667,708,
               715,723,760,763,772,795,928,955,1003,1012,1027,1227,1243,1387,1411,1435,1467,1507,
               1555};
    CNs[8] := {224,260,264,276,308,320,336,360,384,420,456,468,480,504,528,544,552,564,576,580,600,
               612,616,624,640,651,660,672,720,736,768,792,819,820,832,840,852,868,880,900,912,915,
               952,960,987,1008,1032,1035,1060,1092,1120,1128,1131,1152,1155,1204,1240,1248,1275,
               1288,1312,1320,1332,1360,1380,1395,1408,1428,1443,1488,1540,1600,1632,1635,1659,1672,
               1683,1752,1768,1771,1780,1792,1827,1848,1947,1992,1995,2020,2035,2067,2080,2088,2115,
               2128,2139,2163,2212,2272,2275,2368,2392,2451,2475,2632,2667,2715,2755,2788,2832,2907,
               2968,3003,3040,3172,3243,3315,3355,3507,3627,3712,3843,4048,4123,4323,5083,5467,6307};
    CNs[16] := {1056,1140,1344,1440,1560,1680,1716,1824,1860,1872,1920,2016,2040,2100,2112,2176,
                2208,2244,2280,2304,2320,2331,2340,2379,2400,2436,2464,2496,2520,2580,2640,2688,
                2760,2772,2880,3060,3108,3168,3192,3220,3280,3360,3432,3480,3520,3588,3600,3640,
                3648,3795,3808,3828,3840,4020,4032,4128,4180,4260,4275,4368,4420,4440,4452,4480,
                4488,4512,4515,4680,4740,4788,4960,4992,5115,5152,5160,5187,5208,5248,5280,5328,
                5355,5412,5440,5460,5520,5712,5952,6052,6123,6160,6195,6328,6355,6400,6420,6435,
                6528,6580,6612,6688,6708,6820,6840,6867,7008,7035,7072,7120,7315,7392,7395,7480,
                7540,7672,7755,7968,7995,8008,8052,8080,8320,8352,8512,8547,8680,8715,8835,8932,
                9108,9243,9568,9595,9867,9955,10168,10528,10803,10948,11067,11328,11715,11872,
                12160,12483,12595,12688,12915,13195,14008,14155,14547,14763,16192,16555,17427,
                19947,20155};
    CNs[32] := {5760,6240,6360,6720,6900,7140,7488,8160,8400,8448,8580,9120,9240,9280,9540,9600,
                10080,10560,10920,11040,11160,11400,12180,12240,12432,12672,12768,13120,13440,13728,
                13860,13920,14100,14280,14352,14400,14560,14592,14820,15400,16128,16720,17220,17472,
                17760,17952,18720,19152,19240,19320,19380,19635,20020,20148,20475,20608,20640,20832,
                21120,21312,21840,22080,22848,23640,24640,25312,25608,26712,26832,27280,27360,28032,
                28288,28480,28548,29568,29920,30100,30160,30340,30688,31395,32032,32320,33408,33915,
                34720,34840,36432,37107,40672,40755,41475,42112,43435,44115,45312,45747,46852,50752,
                53475,56032,57387,57715,82555};
    CNs[64] := {};
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
    // the congruence alone does not see.  Same test as
    // FieldsOfDefinitionOfCMPointFast, which returns no field in this case.
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
    // m = D_R * N_star_R.  It identifies a point with its conjugate on the
    // quotient exactly when some element of W differs from m by a Galois
    // Atkin-Lehner, which halves the degree.  The product m * mm only has to be
    // Galois; it need not itself lie in W, and requiring that suppresses the
    // halving on quotients such as W = {1, D*N}.
    m := D_R*N_star_R;
    cc_active := exists{mm : mm in W |
        ((D*N) div (D_R*N_R)) mod AtkinLehnerMul(m, mm, D*N) eq 0};

    h_R := #PicardGroup(R);
    // The division is exact: #W_gal / eps_factor is the order of the image of
    // W_gal in Gal(H_R/K), a subgroup of a group of order h_R.  An inexact
    // division means W_gal or the eps correction is wrong, which integer
    // division would otherwise hide by truncating.
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

intrinsic RationalCMDiscs(X::ShimuraQuot) -> Assoc
{Discriminants of the orders whose CM points on X are rational, mapped to the
 number of such points.  Uses Shimura reciprocity through
 DegreeOfFieldOfDefinitionOfCMPoint; requires D*N squarefree.}
    require IsSquarefree(X`D * X`N) : "RationalCMDiscs requires D*N squarefree";
    CNs := CMClassLists();
    cm_pts := AssociativeArray();
    omegaDN := #PrimeFactors(X`D * X`N);
    iterates := [2^i : i in [0..omegaDN]];
    for n in iterates do
        error if not IsDefined(CNs, n),
            Sprintf("CMClassLists() missing class number %o; add CNs[%o] (omega(D*N)=%o, D=%o, N=%o)",
                    n, n, omegaDN, X`D, X`N);
        for a in CNs[n] do
            disc := -a;
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
    CNs := CMClassListsDeg2();
    cm_pts := AssociativeArray();
    omegaDN := #PrimeFactors(X`D * X`N);
    // one power of two further than the rational scan: a degree-2 point can come
    // from an order whose class number is twice as large.
    iterates := [2^i : i in [0..omegaDN+1]];
    if max_class_num gt 0 then
        iterates := [n : n in iterates | n le max_class_num];
    end if;
    for n in iterates do
        error if not IsDefined(CNs, n),
            Sprintf("CMClassListsDeg2() missing class number %o; add CNs[%o] (omega(D*N)=%o, D=%o, N=%o)",
                    n, n, omegaDN, X`D, X`N);
        for a in CNs[n] do
            disc := -a;
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
