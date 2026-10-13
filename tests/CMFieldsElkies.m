// tests/CMFieldsElkies.m
//
// Fields of definition of CM points on X_0^D(1)/W for D = 6, 10, 14, 15 and every W, checked against
// coordinates published independently of this code: Elkies' values of the Hauptmodul t of X*(D) at
// its rational CM points, with the double covers X_0^D(1)/<w_m> -> X*(D) written as y^2 = g_m(t).
// A CM point with t = t0 in Q has field Q(sqrt g_m(t0)) on X/<w_m> (Q when t0 is a branch point),
// the compositum over the three m on X_0^D(1), and Q on X*.  Most of the value is in the rows for
// NON-MAXIMAL orders (conductor noted inline), where the class of complex conjugation in
// Pic(R)/Pic(R)^2 is not determined by the quaternion algebra alone.
//
// SOURCES.
//  [Elk] N. D. Elkies, "Shimura curve computations", ANTS-III, LNCS 1423 (1998); arXiv:math/0005160.
//  [Err] E. Errthum, "Singular moduli of Shimura curves", Canad. J. Math. 63 (2011); arXiv:0711.4316v2.
//  [GR06] J. Gonzalez, V. Rotger, "Non-elliptic Shimura curves of genus one", J. Math. Soc. Japan 58 (2006).
//
//  D = 15: [Elk] Sec. 5.2: the function field of X_0^15(1) is Q(t)(sqrt(-3t), sqrt(-3(t-1)(t-81))).
//          Values and the labels of the elliptic points: [Elk] Table 6 (t = oo, 0, 1, 81 at
//          d = -3, -12, -60, -15).  The prose of Sec. 5.2 lists the discriminants of the elliptic
//          points in a different order from the table; the table is the one consistent with the
//          fixed points below.
//  D = 6:  [Elk] Sec. 3.1: t = 0, 1, oo at the elliptic points of orders 2, 4, 6 (d = -24, -4, -3);
//          the function field of X_0^6(1) is Q(t)(sqrt(-t), sqrt(3(t-1))).  Values: [Elk] Table 1,
//          except d = -40, where Table 1 prints B = 25: [Elk] Table 2 (B = 5^3, |A - B| = 3^7) and
//          [Err] Table "Rational CM points on X_6^*" (s = 5^3) both give t = 2312/125.
//  D = 14: [Elk] Sec. 5.1: t(-4) = oo, t(-8) = 0, the d = -56 pair at 16t^2+13t+8 = 0; the function
//          field of X_0^14(1) is Q(t)(sqrt(-t), sqrt(-(16t^2+13t+8))).  Values: [Elk] Table 5.
//  D = 10: [Elk] Sec. 4, Table 3 (= [Err] Table "Rational CM points on X_10^*", which adds d = -235):
//          t = 0, oo, 2, 27 at d = -3, -8, -20, -40.  [Elk] gives no equations for the covers, so
//          these use the committed models data/models/models_10_1.m pulled back by x = 2/t (the
//          Mobius map matching their branch points to Elkies' elliptic points).  These rows check
//          the field computation given that model; they are not fully independent of this repo.
//
// WHICH COVER IS WHICH.  X/<w_m> -> X* is branched exactly at the images of the fixed points of the
// other two involutions; w_m fixes the CM points by Z[sqrt -m] and, if -m = 1 mod 4, by
// Z[(1+sqrt -m)/2], when these embed in B_D.  D=6: w_2 {-4}, w_3 {-3}, w_6 {-24}.  D=10: w_2 {-8},
// w_5 {-20}, w_10 {-40}.  D=14: w_2 {-4,-8}, w_7 none, w_14 {-56}.  D=15: w_3 {-3,-12}, w_5 none,
// w_15 {-15,-60}.
//
// One point per discriminant: Pic(R) x W acts transitively on the CM points by R of X_0^D(1)
// ([GR06] Prop. 5.6; every conductor below is coprime to D), so X* has exactly one such point,
// the rational one listed.

function cm_norm_val(g, t0)
    QQ := Rationals();
    if t0 cmpeq Infinity() then
        return IsOdd(Degree(g)) select QQ!0 else LeadingCoefficient(g);
    end if;
    k := 0; h := g;
    while Evaluate(h, t0) eq 0 do
        h := h div (Parent(g).1 - t0); k +:= 1;
    end while;
    return IsOdd(k) select QQ!0 else Evaluate(h, t0);
end function;

// Q(sqrt v : v in vals, v ne 0), as an absolute field
function cm_field_from_vals(vals)
    QQ := Rationals();
    sq := [Integers() | ];
    for v in vals do
        if v eq 0 then continue; end if;
        s := SquarefreeFactorization(Numerator(v)*Denominator(v));
        if s ne 1 then Append(~sq, s); end if;
    end for;
    if #sq eq 0 then return QQ; end if;
    K := QuadraticField(sq[1]);
    for s in sq[2..#sq] do
        if not IsSquare(K!s) then
            _<z> := PolynomialRing(K);
            return AbsoluteField(ext<K | z^2 - s>);
        end if;
    end for;
    return K;
end function;

procedure test_CMFieldsElkies()
    QQ := Rationals();
    _<t> := PolynomialRing(QQ);
    // <D, [<m, g_m>], [<d, t0>]>
    data := [*
      // d = -15 comes before every other order of Q(sqrt -15): on Magma V2.29-4, an Artin map
      // built for one order of an imaginary quadratic field can corrupt a later one for another
      // order of the same field (Magma-Maths/Magma#71), and d = -15 then fails after d = -60.
      <15, [<15, -3*t>, <3, -3*(t-1)*(t-81)>, <5, t*(t-1)*(t-81)>],
          [<-15, QQ!81>,                               // elliptic
           <-3, Infinity()>, <-7, QQ!-27>, <-40, 27/2>, <-43, -27/16>, <-195, 81/16>,
           <-120, -243/2>, <-115, 621/16>, <-123, 2025/16>, <-67, -3267/16>,
           <-12, QQ!0>,                                // conductor 2, elliptic
           <-60, QQ!1>,                                // conductor 2, elliptic
           <-48, QQ!243>,                              // conductor 4
           <-28, QQ!-675>,                             // conductor 2
           <-147, -729/112>,                           // conductor 7
           <-240, QQ!9801>]>,                          // conductor 4
      <6, [<2, -t>, <6, 3*(t-1)>, <3, -3*t*(t-1)>],
          [<-3, Infinity()>, <-4, QQ!1>, <-24, QQ!0>, <-84, -169/27>,
           <-40, 2312/125>,                            // [Elk] Table 2 and [Err]; Table 1 misprints
           <-51, -1377/1024>, <-19, 3211/1024>, <-120, 5776/3375>, <-52, 6877/15625>,
           <-132, 13689/15625>, <-168, -701784/15625>, <-43, 21250987/16000000>,
           <-88, 15545888/20796875>, <-312, 27008742384/27680640625>,
           <-75, 152881/138240>,                       // conductor 5
           <-100, 421850521/1771561>,                  // conductor 5
           <-147, -1073152081/3024000000>]>,           // conductor 7
      <14, [<14, -t>, <2, -(16*t^2+13*t+8)>, <7, t*(16*t^2+13*t+8)>],
          [<-4, Infinity()>, <-8, QQ!0>, <-11, QQ!-1>, <-35, QQ!-5>, <-84, 4/9>, <-91, -13/81>,
           <-43, -25/81>, <-168, 200/9>, <-88, -200/81>, <-67, -1225/81>, <-280, -845/1296>,
           <-148, 1225/5184>,
           <-36, 4/3>,                                 // conductor 3
           <-72, 75/16>,                               // conductor 3
           <-100, -196/405>]>,                         // conductor 5
      // covers: data/models/models_10_1.m under x = 2/t (see header)
      <10, [<10, -2*(t-2)>, <5, t-27>, <2, -2*(t-2)*(t-27)>],
          [<-3, QQ!0>, <-8, Infinity()>, <-20, QQ!2>, <-40, QQ!27>, <-52, -54/25>, <-120, -27/49>,
           <-35, 64/7>, <-43, 1728/1225>, <-88, 3375/98>, <-115, 13824/3887>,
           <-280, 35937/7406>, <-67, -216000/8281>,
           <-235, 2^6*3^3*17^3/(7^2*37^2*47)>,         // [Err] only
           <-27, -192/25>,                             // conductor 3
           <-72, 125/147>,                             // conductor 3
           <-180, -2662/169>]>                         // conductor 3
    *];
    printf "Testing CM fields of definition on X_0^D(1)/W against Elkies' coordinates...";
    nchecks := 0;
    for datum in data do
        D, covers, pts := Explode(datum);
        Ws := [{1, c[1]} : c in covers] cat [{1}, Set(Divisors(D))];
        for pt in pts do
            d, t0 := Explode(pt);
            for W in Ws do
                if #W eq 4 then
                    K := QQ;
                elif #W eq 2 then
                    K := cm_field_from_vals([cm_norm_val(c[2], t0) : c in covers | c[1] eq Max(W)]);
                else
                    K := cm_field_from_vals([cm_norm_val(c[2], t0) : c in covers]);
                end if;
                X := CreateShimuraQuot(D, 1, W);
                F := FieldOfDefinitionOfCMPoint(X, d);
                assert IsIsomorphic(AbsoluteField(F), K);
                F := FieldOfDefinitionOfCMPointSlow(X, d);
                assert IsIsomorphic(AbsoluteField(F), K);
                assert DegreeOfFieldOfDefinitionOfCMPoint(X, d) eq Degree(K);
                nchecks +:= 1;
            end for;
        end for;
    end for;
    // every (D, d, W) above, five W per D: (16 + 17 + 15 + 16) * 5
    assert nchecks eq 320;
    printf "Done! (%o checks)\n", nchecks;
end procedure;

test_CMFieldsElkies();
