// The CM values that the special-fibre test uses at D = 10 and D = 22 (special_fiber_cm.m) are those of
// Guo-Yang's Theorem B: the star is X_0^D(1)*, every discriminant is fundamental, and Schofer's formula
// at N = 1 has no m = 0 term.  The test recomputes both tables by that formula (ValuesAtCMPoints) and
// checks, against the stored data:
//   * every stored discriminant is fundamental;
//   * the stored hauptmodul values, rational and of degree two, after the Mobius map fixed by three
//     reference points (the normalisation of a hauptmodul is a choice; every further point is a check);
//   * each cover X_0^D(1)/<w_q> as y^2 = f_q(s): the recomputed y^2 is f_q(s) times a square, so the
//     stored conic is the cover itself and not a quadratic twist of it.
// About a minute.

import "special_fiber_cm.m" : d10cm_discs, d10cm_s, d10cm_fq, d22cm_discs, d22cm_s, d22cm_fq;

procedure check_star_table(DD, discs, svals, fq, curves)
    assert &and[IsFundamentalDiscriminant(d) : d in discs];
    Xstar := rep{X : X in curves | X`D eq DD and X`N eq 1 and IsStarCurve(X)};
    // MaxNum = #discs: the table is exactly the stored discriminants (with the hauptmoduls' zeros and
    // poles, which are among them)
    tab := ValuesAtCMPoints(Xstar, curves : MaxNum := #discs, Prec := 100, Keep := Set(discs));
    assert #tab`AmbiguousSigns eq 0 and IsEmpty(tab`UnscaledKeys);
    ds := tab`Discs; V := tab`Values; si := tab`sIndex;
    assert Set(discs) subset Set(ds);

    Qx<x> := PolynomialRing(Rationals());
    israt := func< v | Type(v) in {FldRatElt, RngIntElt} >;
    // the Mobius map m with m(recomputed) = stored, from three rational reference points
    refs := [d : d in discs | Type(svals[Index(discs, d)]) eq Infty or israt(svals[Index(discs, d)])][1..3];
    pt := func< v | Type(v) eq Infty select [Rationals() | 1, 0] else [Rationals() | v, 1] >;
    R := PolynomialRing(Rationals(), 4);
    eqs := [ (R.1*S[1] + R.2*S[2]) * T[2] - (R.3*S[1] + R.4*S[2]) * T[1]
             where S := pt(V[si][Index(ds, d)]) where T := pt(svals[Index(discs, d)]) : d in refs ];
    K := Kernel(Transpose(Matrix(Rationals(), [[MonomialCoefficient(f, R.j) : j in [1..4]] : f in eqs])));
    assert Dimension(K) eq 1;
    m := Eltseq(Basis(K)[1]);
    mob := func< v | Type(v) eq Infty select (m[3] eq 0 select Infinity() else m[1]/m[3])
                     else (m[3]*v + m[4] eq 0 select Infinity() else (m[1]*v + m[2]) / (m[3]*v + m[4])) >;
    for i -> d in discs do
        if d in refs then continue; end if;
        w := V[si][Index(ds, d)]; st := svals[i];
        if Type(st) eq Infty or israt(st) then
            assert Type(mob(w)) eq Type(st) and (Type(st) eq Infty or mob(w) eq st);
        else
            // a point of degree two, stored as the minimal polynomial of s
            h := Type(w) eq RngUPolElt select Qx!Eltseq(w) else MinimalPolynomial(w);
            num := Numerator(Evaluate(Qx!Eltseq(st), (m[1]*x + m[2]) / (m[3]*x + m[4])));
            assert Degree(h) eq 2 and Normalize(num) eq Normalize(h);
        end if;
    end for;

    for q in Keys(fq) do
        r := [r : r in [1..#tab`Keys_fs] | tab`Keys_fs[r] gt 0 and curves[tab`Keys_fs[r]]`W eq {Integers() | 1, q}][1];
        cf := fq[q]; n := 0;
        for i -> d in discs do
            st := svals[i]; y2 := V[r][Index(ds, d)];
            if not (israt(st) and israt(y2)) or y2 eq 0 then continue; end if;
            fv := cf[1] + cf[2]*st + cf[3]*st^2;
            assert fv ne 0 and IsSquare(y2 / fv);
            n +:= 1;
        end for;
        assert n ge 3;
    end for;
end procedure;

curves := GetHyperellipticCandidates();
check_star_table(10, d10cm_discs, d10cm_s, d10cm_fq, curves);
check_star_table(22, d22cm_discs, d22cm_s, d22cm_fq, curves);
