// FilterByAutomorphismGroup (Brandt-Stichtenoth lemma on the known automorphism group G_Y).
//
// Known values, from the 2026-09-25 prototype (handoff_2026-09-25/scripts/gy.m, reach2.m) and
// reproduced by this implementation over the whole data set (0 of 3423 recorded-subhyperelliptic
// curves of genus >= 2 ruled out; exactly these 9 of the 827 undecided ones):
//   2555 (6,133)/<19,42>, 2568 (6,133)/<21,38>, 4190 (10,27)/<270>, 5635 (14,33)/<2,77>,
//   5639 (14,33)/<6,77>, 6616 (15,28)/<7,60>, 8495 (33,4)/<132>: residual AL group V4, g = 4;
//   7926 (22,45)/<5,9,22>, 7932 (22,45)/<9,10,11>: V4 = <w2, V3>, g = 4.
// 7932 is left out below only for time: deciding that V3 and V3*w2 are not iota takes two
// TraceDNewQuotient calls at level 990 (~30 s each), exactly as for 7926.
// Controls: hyperelliptic curves where some subgroup WOULD fail if iota were not excluded, so the
// iota identification is what saves them, including a non-AL iota (799, iota = V2*w3), and the
// genus-2 curves (1,40)/<w5> and (1,48)/<w3> whose iota is V2*w8, resp. V2*w16 (quotient
// genus 0).

fresh := function(X)
    Y := CreateShimuraQuot(X`D, X`N, X`W);
    Y`g := X`g;
    Y`CurveID := X`CurveID;
    return Y;
end function;

procedure test_AutomorphismGroupRuledOut(curves)
    known := [<2555, "H = <w2, w3> is V4 and g = 4 is even">,
              <2568, "H = <w2, w3> is V4 and g = 4 is even">,
              <4190, "H = <w2, w10> is V4 and g = 4 is even">,
              <5635, "H = <w7, w3> is V4 and g = 4 is even">,
              <5639, "H = <w7, w2> is V4 and g = 4 is even">,
              <6616, "H = <w3, w4> is V4 and g = 4 is even">,
              <8495, "H = <w11, w3> is V4 and g = 4 is even">,
              <7926, "H = <w2, V3> is V4 and g = 4 is even">];
    cs := [fresh(curves[k[1]]) : k in known];
    FilterByAutomorphismGroup(~cs);
    for i->k in known do
        X := cs[i];
        assert assigned X`IsSubhyp and not X`IsSubhyp and not X`IsHyp;
        printf "  %o: %o\n", X`CurveID, X`TestInWhichProved;
        assert X`TestInWhichProved eq "AutomorphismGroup " cat k[2];
    end for;
end procedure;

procedure test_AutomorphismGroupControls(curves)
    // iota = V2*W_(2^a) on the genus-2 curves (1,40)/<w5> (149, 2^a = 8) and (1,48)/<w3>
    // (191, 2^a = 16), as predicted by the lemma for G_Y = D4 = <S2, W_(2^a)>
    for c in [<149, 40, "V2*w8">, <191, 48, "V2*w16">] do
        X := curves[c[1]];
        assert <X`D, X`N, X`g> eq <1, c[2], 2>;
        G, _ := KnownAutomorphismGroup(X);
        assert IdentifyGroup(G) eq <8, 3>;
        iotas := PossibleHyperellipticInvolutions(X);
        printf "  %o: possible iota %o\n", c[1], iotas;
        assert iotas eq [c[3]];
        assert CheckAutomorphismGroup(X);
    end for;
    // recorded hyperelliptic, genus >= 3, and some subgroup contains iota:
    //   242 (1,60)/<12> g=4 iota w15, 348 (1,78)/<6> g=6 iota w26, 799 (1,168)/<21,56> g=4 iota V2*w3,
    //   322 (1,72)/<9> g=3 iota V2*V3*w8, 146 (1,40) g=3 iota V2*w40, 613 (1,126)/<7,9> g=3 iota V3
    ctl := [<242, "w15">, <348, "w26">, <799, "V2*w3">, <322, "V2*V3*w8">, <146, "V2*w40">, <613, "V3">];
    for c in ctl do
        assert curves[c[1]]`IsSubhyp;
        assert c[2] in PossibleHyperellipticInvolutions(curves[c[1]]);
    end for;
    cs := [fresh(curves[c[1]]) : c in ctl];
    FilterByAutomorphismGroup(~cs);
    assert &and[not assigned X`IsSubhyp : X in cs];
    // decided curves are skipped (ShimuraQuot is a reference type, so compare against a saved string)
    proof := curves[242]`TestInWhichProved;
    cs := [curves[242]];
    FilterByAutomorphismGroup(~cs);
    assert cs[1]`IsSubhyp and cs[1]`TestInWhichProved eq proof;
end procedure;

curves := GetHyperellipticCandidates();
test_AutomorphismGroupRuledOut(curves);
test_AutomorphismGroupControls(curves);
