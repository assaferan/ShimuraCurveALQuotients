// FilterByTwistedWeilPolynomial: the Weil polynomial P_+(t) P_-(-t) of the twist by an involution h
// defined over Q must be in the LMFDB hyperelliptic tables data/hypg<g>q<p>.txt.
//
// Known values, from the twisted-Weil sweep (sweeps/twisted_weil/tw.m on branch twisted-trace-sweep;
// the AL hits were cross-checked by Eichler-Selberg):
//   1583 X_0^6(23)/<w23>   fails at p = 5 for h = w2 and h = w6
//   420  X_0(90)/<w2,w9>   fails at p = 13 for h = w5
// At h = 1 the code must reproduce the pipeline's own Weil polynomials: on 1583 it equals
// WeilPolynomial (trace formula, no modular symbols), and on 169 (a recorded "WeilPolynomial with
// p = 2" verdict) it is not in the table.  (The sweep reproduced all 60 recorded WeilPolynomial
// verdicts this way.)  Controls: recorded-hyperelliptic curves with AL and V3 twists.

fresh := function(X)
    Y := CreateShimuraQuot(X`D, X`N, X`W);
    Y`g := X`g;
    Y`CurveID := X`CurveID;
    return Y;
end function;

procedure test_TwistedWeil(curves)
    X := curves[1583];
    assert <X`D, X`N, X`g> eq <6, 23, 3>;
    tw := TwistedWeilPolynomials(X, 5);
    assert tw[1][1] eq "1" and tw[1][2] eq WeilPolynomial(X, 5);
    _<t> := Universe([x[2] : x in tw]);
    P := AssociativeArray(); for x in tw do P[x[1]] := x[2]; end for;
    assert P["w2"] eq t^6 - 4*t^5 + 19*t^4 - 40*t^3 + 95*t^2 - 100*t + 125;
    assert P["w6"] eq Evaluate(P["w2"], -t);
    // h = 1 reproduces a recorded WeilPolynomial verdict
    Y := curves[169];
    assert Y`TestInWhichProved eq Sprintf("WeilPolynomial with p = %o", 2);
    P1 := TwistedWeilPolynomials(Y, 2)[1][2];
    tab := {l : l in Split(Read(Sprintf("data/hypg%oq%o.txt", Y`g, 2)), "\n") | #l gt 0};
    assert "[" cat Join([IntegerToString(x) : x in Reverse(Coefficients(P1))], ",") cat "]" notin tab;
    // the filter
    known := [<1583, "w2", 5>, <420, "w5", 13>];
    ctl := [1577, 1579, 263, 613];
    assert &and[curves[c]`IsSubhyp : c in ctl];
    cs := [fresh(curves[k[1]]) : k in known] cat [fresh(curves[c]) : c in ctl];
    FilterByTwistedWeilPolynomial(~cs);
    for i->k in known do
        X := cs[i];
        printf "  %o: %o\n", X`CurveID, assigned X`TestInWhichProved select X`TestInWhichProved else "-";
        assert assigned X`IsSubhyp and not X`IsSubhyp and not X`IsHyp;
        assert X`TestInWhichProved eq Sprintf("TwistedWeilPolynomial h = %o with p = %o", k[2], k[3]);
    end for;
    assert &and[not assigned cs[i]`IsSubhyp : i in [#known+1..#cs]];
end procedure;

// Star curves (FilterByTwistedWeilPolynomialStar): X_0^(39)(16)^* (9698, undecided) fails at p = 5
// for h = V2; the hyperelliptic star curves X_0^(35)(16)^* (9244) and X_0^*(176) (828), both with
// V2, do not.
procedure test_TwistedWeilStar(curves)
    X := curves[9698];
    assert <X`D, X`N> eq <39, 16> and IsStarCurve(X);
    ctl := [9244, 828];
    assert &and[IsStarCurve(curves[c]) and curves[c]`IsSubhyp : c in ctl];
    cs := [fresh(X)] cat [fresh(curves[c]) : c in ctl];
    FilterByTwistedWeilPolynomial(~cs);
    printf "  star %o: %o\n", X`CurveID, cs[1]`TestInWhichProved;
    assert (not cs[1]`IsSubhyp) and cs[1]`TestInWhichProved eq "TwistedWeilPolynomial h = V2 with p = 5";
    assert &and[not assigned cs[i]`IsSubhyp : i in [2..#cs]];
end procedure;

curves := GetHyperellipticCandidates();
test_TwistedWeil(curves);
test_TwistedWeilStar(curves);
