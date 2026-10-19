// The positions of P_-4 and P_-8 on the genus-one stars X_0^142(1)* and X_0^158(1)*, checked
// WITHOUT any Borcherds form or Schofer value.
//
// Those two positions are what fixes the branch loci recorded in data/models/PROVENANCE.md, and
// at 158 the model is the only source for them, so an independent check is worth having.
//
// The argument.  w_2 fixes the CM points of discriminants -4 and -8 and w_{D} fixes none, so
// C = X/<w_{2D}> -> X* is a double cover of the star branched exactly over P_-4 and P_-8; by
// Riemann-Hurwitz C has genus 2, and Jac(C) = Jac(X*) x Prym with Prym an elliptic curve.  By
// Jacquet-Langlands the quaternionic eigenvalue of w_p on a newform is its p-th coefficient, so
// the star is the rational newform with Atkin-Lehner signs (-1, -1) -- that is, 142a1 and 158a1 --
// and the other factor of Jac(C), the one with quaternionic w_{2D}-eigenvalue +1 that is not the
// star, is the rational newform with Atkin-Lehner signs (+1, +1).  The Prym is that newform up to
// a quadratic twist ramified only where C and the newform can differ, i.e. at the primes of 2D.
//
// Translating the branch pair by one of its points, C is the double cover of E branched over O
// and n*Q, where n*Q = P_-4 - P_-8 and Q generates E(Q) = Z; translation is an automorphism of E,
// so this does not change C.  Such a cover is y^2 = h with div h = (O) + (nQ) - 2(P) and
// 2P = nQ, so it is defined over Q only for even n -- as it must be, C being defined over Q.
// Running over the even n, the Prym matches the (+1, +1) newform for exactly one of them, and it
// is the n the models record: 2 at 142 (where Guo-Yang's table gives P_-4 = -2Q, P_-8 = O) and 6
// at 158 (where the model gives P_-4 = O, P_-8 = 6Q).
//
// Sources: the newforms and their Atkin-Lehner signs are computed here from the modular symbols
// of each level; the curves 142a1 and 158a1 and the generator Q are the stored models; the
// positions n = 2 and n = 6 are what is being checked.

procedure test_prym_branch_positions()
    printf "Branch positions on the genus-one stars from Prym traces (no Borcherds form)...";
    // base: level, the stored model, the generator of E(Q), and |P_-4 - P_-8| in units of Q
    cases := [<142, [1,-1,1,-12,15], [1,1,1], 2>, <158, [1,-1,1,-9,9], [-1,-4,1], 6>];
    for cs in cases do
        D, aInv, Qc, expected := Explode(cs);
        E := EllipticCurve([Rationals() | a : a in aInv]);
        Q := E!Qc;
        assert #TorsionSubgroup(E) eq 1 and Rank(E) eq 1;
        // the rational newform of level D with Atkin-Lehner signs (+1, +1); at q || D the sign of
        // a rational newform is -a_q
        target := [* *];
        for orb in Newforms(CuspForms(Gamma0(D), 2)) do
            f := orb[1];
            if Degree(BaseRing(Parent(f))) ne 1 then continue; end if;
            if forall{q : q in PrimeDivisors(D) | -Coefficient(f, q) eq 1} then Append(~target, f); end if;
        end for;
        assert #target eq 1;
        f := target[1];
        twists := [c : c in [-2*D..2*D] | c ne 0 and IsSquarefree(c)
                                         and forall{q : q in PrimeDivisors(Abs(c)) | D mod q eq 0 or q eq 2}];
        FE := FunctionField(E);
        good := [p : p in PrimesUpTo(47) | D mod p ne 0];
        matched := [];
        for n in [2, 4, 6, 8, 10] do
            P := (n div 2)*Q;
            V, mV := RiemannRochSpace(2*Place(P) - Place(E!0) - Place(n*Q));
            assert Dimension(V) eq 1;                  // the cover exists for every even n
            h := mV(V.1);
            pf := ProjectiveFunction(h);
            Nh := Numerator(pf); Dh := Denominator(pf);
            den := LCM([Denominator(x) : x in Coefficients(Nh) cat Coefficients(Dh)]);
            aP := AssociativeArray();
            for p in good do
                Ep := ChangeRing(E, GF(p));
                FEp := FunctionField(Ep);
                AFFp, toAFFp := AlgorithmicFunctionField(FEp);
                R3 := PolynomialRing(GF(p), 3);
                xp := FEp.1; yp := FEp.2;
                Np := &+[R3 | (GF(p)!(den*MonomialCoefficient(Nh, m)))*Monomial(R3, Exponents(m)) : m in Monomials(Nh)];
                Dp := &+[R3 | (GF(p)!(den*MonomialCoefficient(Dh, m)))*Monomial(R3, Exponents(m)) : m in Monomials(Dh)];
                if Evaluate(Dp, [xp, yp, 1]) eq 0 or Evaluate(Np, [xp, yp, 1]) eq 0 then continue; end if;
                hp := Evaluate(Np, [xp, yp, 1]) / Evaluate(Dp, [xp, yp, 1]);
                PTp<Tp> := PolynomialRing(AFFp);
                try
                    FC := FunctionField(Tp^2 - toAFFp(hp));
                    if Genus(FC) ne 2 then continue; end if;
                    aC := p + 1 - #Places(FC, 1);
                catch e
                    continue;
                end try;
                aP[p] := aC - TraceOfFrobenius(Ep);     // Jac(C) = E x Prym
            end for;
            used := Sort([p : p in Keys(aP)]);
            assert #used ge 12;
            // the scaling of h is arbitrary up to a square class, so the twist is quantified over
            if exists{c : c in twists | forall{p : p in used | KroneckerSymbol(c, p)*aP[p] eq Coefficient(f, p)}} then
                Append(~matched, n);
            end if;
        end for;
        assert matched eq [expected];
    end for;
    printf "Done!\n";
end procedure;

test_prym_branch_positions();
