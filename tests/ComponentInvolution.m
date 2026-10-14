// Criterion H (ComponentInvolution.m): at an odd p || N with p !| D, a special-fibre component of
// genus >= 2 whose hyperelliptic involution is an Atkin-Lehner w_u forces iota = w_u.
// (a) The 12 open curves of the 2026-09-30 list that criterion H rules out, from the prototype run.
//     All 12 also fail an independent Brandt-module check of the node gluing under w_u, and the
//     first four are also ruled out by a point-count criterion at p (both computed separately).
// (b) Soundness: no curve proved subhyperelliptic by other stages (genus >= 2) is ruled out.
// (c) A genus-1 component must not be used.  (d) In case B the component is X_0^D(N/p)/W',
//     W' = {w in W : p !| w}, not the image of W.

procedure test_ComponentInvolutionRules()
    printf "Testing ComponentInvolutionNotHyperelliptic on the twelve open curves...";
    cases := [ <14, 15, {1,3,7,21}>, <14, 15, {1,3,70,210}>, <15, 28, {1,7,12,84}>, <15, 28, {1,3,7,21}>,
               <14, 15, {1,5,7,35}>, <15, 14, {1,2,7,14}>, <21, 20, {1,15,28,420}>,
               <21, 20, {1,15,84,140}>, <21, 20, {1,12,15,20,21,28,35,420}>,
               <22, 15, {1,2,55,110}>, <34, 21, {1,2,7,14,51,102,357,714}>, <62, 5, {1,5,31,155}> ];
    for c in cases do
        ok, wit := ComponentInvolutionNotHyperelliptic(c[1], c[2], {Integers() | x : x in c[3]});
        printf "\n  D=%o N=%o W=%o: %o %o", c[1], c[2], c[3], ok, ok select wit else "";
        assert ok;
    end for;
    printf "\nDone\n";
end procedure;

procedure test_ComponentInvolutionSound()
    printf "Testing ComponentInvolutionNotHyperelliptic against known subhyperelliptic curves...";
    reached := 0;
    for X in GetHyperellipticCandidates() do
        if not (assigned X`IsSubhyp and X`IsSubhyp) or X`g lt 2 then continue; end if;
        reached +:= 1;
        ok, wit := ComponentInvolutionNotHyperelliptic(X`D, X`N, X`W);
        if ok then printf "\n  *** D=%o N=%o W=%o known subhyperelliptic but ruled out: %o", X`D, X`N, X`W, wit; end if;
        assert not ok;
    end for;
    printf " %o curves. Done\n", reached;
    assert reached gt 0;
end procedure;

procedure test_ComponentInvolutionHypotheses()
    printf "Testing that genus-1 components and the image group are not used...";
    star := func<a, b | a*b div GCD(a, b)^2>;
    // Each curve is hyperelliptic because X/w_v has genus 0 (checked here by the genus formula).
    // (c) X_0(60)/w_60, genus 3, v = 4.  At p = 3, W' = {1}: the component X_0(20) has genus 1,
    //     X_0(20)/w_20 genus 0 and g(X/w_20) = 1, so allowing genus 1 would rule this curve out.
    W := {Integers() | 1, 60};
    assert GenusShimuraCurveQuotient(1, 60, W) eq 3;
    assert GenusShimuraCurveQuotient(1, 60, W join {star(4, w) : w in W}) eq 0;
    assert GenusShimuraCurveQuotient(1, 20, {Integers() | 1}) eq 1;
    assert GenusShimuraCurveQuotient(1, 20, {Integers() | 1, 20}) eq 0;
    assert GenusShimuraCurveQuotient(1, 60, W join {star(20, w) : w in W}) eq 1;
    assert not ComponentInvolutionNotHyperelliptic(1, 60, W);
    // (d) X_0^10(33)/<w_5, w_33>, genus 3, v = 6.  At p = 3, W' = {1, 5}; the image of W is
    //     {1, 5, 11, 55}, whose component has genus >= 2 and a w_2 with genus-0 quotient, and
    //     g(X/w_2) = 1, so using the image would rule this curve out.
    W := {Integers() | 1, 5, 33, 165};
    assert GenusShimuraCurveQuotient(10, 33, W) eq 3;
    assert GenusShimuraCurveQuotient(10, 33, W join {star(6, w) : w in W}) eq 0;
    Wi := {Integers() | 1, 5, 11, 55};
    assert GenusShimuraCurveQuotient(10, 11, Wi) ge 2;
    assert GenusShimuraCurveQuotient(10, 11, Wi join {star(2, w) : w in Wi}) eq 0;
    assert GenusShimuraCurveQuotient(10, 33, W join {star(2, w) : w in W}) eq 1;
    assert not ComponentInvolutionNotHyperelliptic(10, 33, W);
    printf "Done\n";
end procedure;

test_ComponentInvolutionRules();
test_ComponentInvolutionSound();
test_ComponentInvolutionHypotheses();
