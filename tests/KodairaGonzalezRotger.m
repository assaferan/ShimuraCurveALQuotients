// The Kodaira index from the metric dual graph (KodairaIndexFromDualGraph) against the identifications
// of Gonzalez-Rotger, Non-elliptic Shimura curves of genus one, J. Math. Soc. Japan 58 (2006),
// arXiv:math/0612732: the genus-one curves X_0^D(1) of the table in the proof of Theorem 3.3 and the
// genus-one quotients X_0^D(1)/<w_m> of Table 2 of Theorem 4.2. Theirs is the method this function
// implements, so agreement on their curves is the outside check of it.
//
// Each row is tested in two independent steps, and the isogeny class is never taken from a label:
//   [1] the class comes from the trace formula. Jac(X_0^D(1)/W) is D-new, so its conductor is exactly D
//       (Gonzalez-Rotger, Theorem 5.4 and Section 4), and ComputePointsViaTrace gives its traces of
//       Frobenius at the good primes; exactly one isogeny class of conductor D must match all of them.
//   [2] the indices from the dual graph at the primes of D must single out ONE curve of that class.
// The published labels are then a cross-check on [2] where the repository has transcribed them:
// the six full curves and the seventeen quotients of Table 2 are in tests/GonzalezRotger.m (PART 1,
// GR, and PART 4, GR2), and three more quotients are in its PART 2b companion rows. The three
// remaining quotients of p.8 (D = 21, 33, 34 by w_3, w_11, w_2) have no label transcribed in this
// repository, so for them only [1] and [2] run -- which is the whole content of the check anyway.

ROWS := [* // <D, m (1 = the full curve), published Cremona label or "">
    <14, 1, "14a2">, <15, 1, "15a1">, <21, 1, "21a2">,      // GR table in the proof of Theorem 3.3
    <33, 1, "33a1">, <34, 1, "34a3">, <46, 1, "46a2">,
    <14, 7, "14a1">, <15, 5, "15a2">, <46, 23, "46a1">,     // GR p.8, transcribed in PART 2b
    <21, 3, "">, <33, 11, "">, <34, 2, "">,                 // GR p.8, label not transcribed here
    <39, 13, "39a1">, <55, 5, "55a1">, <62, 2, "62a3">,     // GR Table 2
    <69, 3, "69a2">, <77, 11, "77c2">, <85, 17, "85a1">,
    <94, 2, "94a2">, <178, 89, "178b1">, <210, 30, "210e5">,
    <210, 42, "210a6">, <210, 70, "210c3">, <210, 105, "210b6">,
    <330, 3, "330b3">, <330, 22, "330c3">, <330, 33, "330d2">,
    <330, 165, "330a2">, <462, 154, "462b2"> *];

error if #ROWS ne 29,
    Sprintf("the calibration lists %o curves, not the 29 of Gonzalez-Rotger's two tables", #ROWS);

DB := CremonaDatabase();
kgr_bad := [];
kgr_n := 0;
for kgr_t in ROWS do
    D := kgr_t[1]; m := kgr_t[2]; lab := kgr_t[3];
    W := m eq 1 select {Integers() | 1} else {Integers() | 1, m};
    X := CreateShimuraQuot(D, 1, W);
    X`g := GenusShimuraCurveQuotient(D, 1, W);
    error if X`g ne 1, Sprintf("%o/%o has genus %o, not one", D, W, X`g);
    // [1] the isogeny class from the trace formula
    ps := [p : p in PrimesInInterval(3, 60) | D mod p ne 0];
    aps := [p + 1 - ComputePointsViaTrace(X, p, 1) : p in ps];
    cls := [];
    for i in [1..NumberOfIsogenyClasses(DB, D)] do
        Ei := EllipticCurve(DB, D, i, 1);
        if &and[TraceOfFrobenius(Ei, ps[k]) eq aps[k] : k in [1..#ps]] then
            cls := [EllipticCurve(DB, D, i, j) : j in [1..NumberOfCurves(DB, D, i)]];
        end if;
    end for;
    if #cls eq 0 then
        Append(~kgr_bad, Sprintf("%o/%o: no isogeny class of conductor %o has the traces of the "
                                 * "trace formula", D, Sort(SetToSequence(W)), D));
        continue;
    end if;
    // [2] the indices from the dual graph, and which curves of the class have them
    ns := [];
    for p in PrimeDivisors(D) do
        Append(~ns, <p, KodairaIndexFromDualGraph(D, 1, W, p)>);
    end for;
    fits := [CremonaReference(Ec) : Ec in cls |
             &and[KodairaSymbol(Ec, t[1]) eq KodairaSymbol("I" cat IntegerToString(t[2])) : t in ns]];
    kgr_n +:= 1;
    if #fits ne 1 then
        Append(~kgr_bad, Sprintf("%o/%o: indices %o fit %o curves of the class %o, not one",
                                 D, Sort(SetToSequence(W)), ns, #fits, [CremonaReference(Ec) : Ec in cls]));
    elif lab ne "" and fits[1] ne lab then
        Append(~kgr_bad, Sprintf("%o/%o: indices %o give %o, the paper gives %o",
                                 D, Sort(SetToSequence(W)), ns, fits[1], lab));
    end if;
end for;

error if kgr_n ne 29, Sprintf("only %o of the 29 curves were reached", kgr_n);
error if not IsEmpty(kgr_bad),
    Sprintf("Kodaira index against Gonzalez-Rotger: %o failure(s): %o", #kgr_bad, kgr_bad);
printf "Kodaira index from the dual graph against Gonzalez-Rotger: %o curves, each identified "
    * "uniquely in its isogeny class, and each agreeing with the published label where this "
    * "repository has one\n", kgr_n;
