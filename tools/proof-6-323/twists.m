// The class 1938h of the star of X_0^6(323): torsion over Q and over K = Q(sqrt -2) with its minus part, component groups,
// whether the 2-torsion point lies on the identity component at 2 and 3, and the twists of E by -8 and 969.
//   magma -b tools/proof-6-323/twists.m
// a_p of the star from the trace formula (tools/proof-6-323/isogenyclass.m checks them to the Sturm bound)
aps := [ <5, -2>, <7, -2>, <11, -4>, <13, 0>, <23, 0>, <29, 4>, <31, -4>, <37, -10>, <41, -2>, <43, -12>, <47, -2>, <53, -6>, <59, 12> ];
K<s> := QuadraticField(-2);
for E in IsogenousCurves(EllipticCurve("1938h1")) do
    lab := CremonaReference(E);
    ok := forall{t : t in aps | TraceOfFrobenius(E, t[1]) eq t[2]};
    EK := BaseChange(E, K);
    T, m := TorsionSubgroup(E); TK, mK := TorsionSubgroup(EK);
    minus := [mK(g) : g in TK | mK(g) ne EK!0 and EK![Conjugate(c) : c in Eltseq(mK(g))] eq -mK(g)];
    printf "%o: a_p match %o | %o | tors/Q %o | tors/K %o, minus part %o | Kodaira at 2,3,17,19 %o | Tamagawa %o\n",
        lab, ok, aInvariants(E), Invariants(T), Invariants(TK), minus,
        [KodairaSymbol(E, p) : p in [2,3,17,19]], [TamagawaNumber(E, p) : p in [2,3,17,19]];
end for;
E := EllipticCurve("1938h1");
for d in [-8, 969] do
    Ed := MinimalModel(QuadraticTwist(E, d));
    r, proved := Rank(Ed);
    printf "twist by %o: root number %o, rank %o (proved %o), generators %o\n", d, RootNumber(Ed), r, proved, Generators(Ed);
end for;
function InE0(E, P, p)   // P reduces to the identity component: to O when p divides a denominator, else to a smooth point
    if exists{c : c in Eltseq(P) | Valuation(Denominator(c), p) gt 0} then return true; end if;
    Ep := ChangeRing(E, GF(p)); sing := SingularPoints(Curve(Ep));
    return not (Curve(Ep)![GF(p)!c : c in Eltseq(P)]) in sing;
end function;
for lab in ["1938h1", "1938h2"] do
    E := MinimalModel(EllipticCurve(lab)); G, m := TorsionSubgroup(E);
    for g in G do if g ne G!0 then T := m(g);
        printf "%o: T = %o; reduction type at 2, 3: %o; T in E^0 at 2: %o, at 3: %o\n", lab, T, [ReductionType(E, p) : p in [2, 3]], InE0(E, T, 2), InE0(E, T, 3);
    end if; end for;
end for;
quit;
