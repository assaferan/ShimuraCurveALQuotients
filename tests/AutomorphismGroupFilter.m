// FilterByAutomorphismGroup (Brandt-Stichtenoth lemma on the known automorphism group G_Y).
//
// Known values, from the 2026-09-25 prototype (handoff_2026-09-25/scripts/gy.m, reach2.m) and
// reproduced by this implementation over the whole data set (0 of 3423 recorded-subhyperelliptic
// curves of genus >= 2 ruled out; exactly these 9 of the 827 undecided ones):
//   (6,133)/<19,42>, (6,133)/<21,38>, (10,27)/<270>, (14,33)/<2,77>,
//   (14,33)/<6,77>, (15,28)/<7,60>, (33,4)/<132>: residual AL group V4, g = 4;
//   (22,45)/<5,9,22>, (22,45)/<9,10,11>: V4 = <w2, V3>, g = 4.
// (22,45)/<9,10,11> is left out below only for time: deciding that V3 and V3*w2 are not iota takes
// two TraceDNewQuotient calls at level 990 (~30 s each), exactly as for (22,45)/<5,9,22>.
// Controls: hyperelliptic curves where some subgroup WOULD fail if iota were not excluded, so the
// iota identification is what saves them, including a non-AL iota ((1,168)/<21,56>, iota = V2*w3), and the
// genus-2 curves (1,40)/<w5> and (1,48)/<w3> whose iota is V2*w8, resp. V2*w16 (quotient
// genus 0).

// Each curve is named by D, N and its Atkin-Lehner group W: its position in the list of quotients
// changes whenever the enumeration of pairs does.
curve := function(curves, D, N, W)
    found := [X : X in curves | X`D eq D and X`N eq N and X`W eq W];
    assert #found eq 1;
    return found[1];
end function;

fresh := function(X)
    Y := CreateShimuraQuot(X`D, X`N, X`W);
    Y`g := X`g;
    Y`CurveID := X`CurveID;
    return Y;
end function;

procedure test_AutomorphismGroupRuledOut(curves)
    known := [<<6, 133, {1, 19, 42, 798}>, "H = <w2, w3> is V4 and g = 4 is even">,
              <<6, 133, {1, 21, 38, 798}>, "H = <w2, w3> is V4 and g = 4 is even">,
              <<10, 27, {1, 270}>, "H = <w2, w10> is V4 and g = 4 is even">,
              <<14, 33, {1, 2, 77, 154}>, "H = <w7, w3> is V4 and g = 4 is even">,
              <<14, 33, {1, 6, 77, 462}>, "H = <w7, w2> is V4 and g = 4 is even">,
              <<15, 28, {1, 7, 60, 420}>, "H = <w3, w4> is V4 and g = 4 is even">,
              <<33, 4, {1, 132}>, "H = <w11, w3> is V4 and g = 4 is even">,
              <<22, 45, {1, 5, 9, 22, 45, 110, 198, 990}>, "H = <w2, V3> is V4 and g = 4 is even">];
    cs := [fresh(curve(curves, k[1][1], k[1][2], k[1][3])) : k in known];
    FilterByAutomorphismGroup(~cs);
    for i->k in known do
        X := cs[i];
        assert assigned X`IsSubhyp and not X`IsSubhyp and not X`IsHyp;
        printf "  %o: %o\n", X`CurveID, X`TestInWhichProved;
        assert X`TestInWhichProved eq "AutomorphismGroup " cat k[2];
    end for;
end procedure;

procedure test_AutomorphismGroupControls(curves)
    // iota = V2*W_(2^a) on the genus-2 curves (1,40)/<w5> (2^a = 8) and (1,48)/<w3>
    // (2^a = 16), as predicted by the lemma for G_Y = D4 = <S2, W_(2^a)>
    for c in [<{1, 5}, 40, "V2*w8">, <{1, 3}, 48, "V2*w16">] do
        X := curve(curves, 1, c[2], c[1]);
        assert X`g eq 2;
        G, _ := KnownAutomorphismGroup(X);
        assert IdentifyGroup(G) eq <8, 3>;
        iotas := PossibleHyperellipticInvolutions(X);
        printf "  X_0(%o)/%o: possible iota %o\n", c[2], c[1], iotas;
        assert iotas eq [c[3]];
        assert CheckAutomorphismGroup(X);
    end for;
    // recorded hyperelliptic, genus >= 3, and some subgroup contains iota:
    //   (1,60)/<12> g=4 iota w15, (1,78)/<6> g=6 iota w26, (1,168)/<21,56> g=4 iota V2*w3,
    //   (1,72)/<9> g=3 iota V2*V3*w8, (1,40) g=3 iota V2*w40, (1,126)/<7,9> g=3 iota V3
    ctl := [<curve(curves, 1, 60, {1, 12}), "w15">, <curve(curves, 1, 78, {1, 6}), "w26">,
            <curve(curves, 1, 168, {1, 21, 24, 56}), "V2*w3">, <curve(curves, 1, 72, {1, 9}), "V2*V3*w8">,
            <curve(curves, 1, 40, {1}), "V2*w40">, <curve(curves, 1, 126, {1, 7, 9, 63}), "V3">];
    for t in ctl do
        assert t[1]`IsSubhyp;
        assert t[2] in PossibleHyperellipticInvolutions(t[1]);
    end for;
    cs := [fresh(t[1]) : t in ctl];
    FilterByAutomorphismGroup(~cs);
    assert &and[not assigned X`IsSubhyp : X in cs];
    // decided curves are skipped (ShimuraQuot is a reference type, so compare against a saved string)
    proof := ctl[1][1]`TestInWhichProved;
    cs := [ctl[1][1]];
    FilterByAutomorphismGroup(~cs);
    assert cs[1]`IsSubhyp and cs[1]`TestInWhichProved eq proof;
end procedure;

// The witness does not depend on the Magma process.  Before the subgroups were ordered by their
// elements, five runs gave two different witnesses on X_0(72) and four on X_0(120)/<w15>
// (different subgroups of order 8 and different generators).  The strings are a regression pin of this
// implementation; the claim behind each, (Z/2)^3 not in PGL_2, is the lemma.
procedure test_AutomorphismGroupWitnessIsCanonical()
    for c in [<72, {1}, "H = <w8, w9, V2> is <8, 5>, not in PGL_2">,
              <120, {1, 15}, "H = <w3, w8, V2> is <8, 5>, not in PGL_2">] do
        X := CreateShimuraQuot(1, c[1], c[2]);
        X`g := GenusShimuraCurveQuotient(1, c[1], c[2]);
        assert X`g eq 5;
        for i in [1..3] do
            ok, why := CheckAutomorphismGroup(X);
            assert not ok and why eq c[3];
        end for;
        printf "  X_0(%o)/%o: %o\n", c[1], c[2], why;
    end for;
end procedure;

curves := GetHyperellipticCandidates();
test_AutomorphismGroupRuledOut(curves);
test_AutomorphismGroupControls(curves);
test_AutomorphismGroupWitnessIsCanonical();
