// tests/QuotientByInvolution.m -- the helper in tests/_quotbyinvol.m, checked forwards and
// BACKWARDS: every derivation it gets right, and every wrong input it must refuse.
//
// WHY THE CONTROLS COME FIRST.  This helper's output is destined to become the expected value of
// other tests, so an instrument that cannot go red would launder a wrong curve into forty of them.
// PART 1 therefore feeds it five inputs that must each be REJECTED, and the count guard at the
// bottom fails if any control stops firing.  "A passing check is not evidence" -- see
// tests/ConicClasses.m for the same pattern.
//
// SELF-CONTAINED BY DESIGN.  The two bases below carry PUBLISHED equations and PUBLISHED
// involutions (Gonzalez-Rotger, "Non-elliptic Shimura curves of genus one", J. Math. Soc. Japan 58
// (2006) 927-948, Table 1 p.8 and the involution table under it), transcribed here rather than
// read out of tests/X0_14_1.m and tests/X0_6_7.m, so this test does not depend on those files'
// internals.  The cross-base sweep over all 50 X0_ bases is tests/_quotsweep.m, hand-run.
//
// ⚠ COMPARISON IS BY ISOMORPHISM CLASS, NEVER COEFFICIENT-WISE.  The derived model is isomorphic
// to the committed one, not equal: at 14_1, w_2 yields -u^2-13u-128 against a committed
// -u^2+13u-128, because u = x*sigma(x) = -x^2 differs from the committed normalisation by u -> -u.

import "tests/_quotbyinvol.m" : QuotientByInvolution, QuotientMatches,
                                ALMatrixGroupFromGenerators, ALCompose;

Q := Rationals();
P<x> := PolynomialRing(Q);
nControl := 0;   // negative controls that fired
nDeriv   := 0;   // quotients derived and matched

// =============================================================== PART 1: NEGATIVE CONTROLS
// X_0(14,1): y^2 = -x^4+13x^2-128, w_2 = (-x,y), w_7 = (-x,-y), w_14 = (x,-y).
f14 := -x^4 + 13*x^2 - 128;

// The baseline must SUCCEED, or the controls below prove nothing (a broken helper refuses
// everything and would look like five passing controls).
ok, _, _ := QuotientByInvolution(f14, DiagonalMatrix(Q, [-1,1,1]));
error if not ok, "QuotientByInvolution: baseline w_2 at 14_1 failed; controls below are vacuous";

controls := [* *];
//  (1) an off-diagonal perturbation: sigma is still an involution, but the map no longer
//      preserves the curve.  This is the control that matters -- a wrong involution otherwise
//      yields a perfectly plausible curve.
Append(~controls, <"perturbed matrix", Matrix(Q, 3,3, [-1,0,0,  0,1,0,  1,0,1])>);
//  (2) not an involution at all
Append(~controls, <"non-involution",   Matrix(Q, 3,3, [ 2,0,0,  0,1,0,  0,0,1])>);
//  (3) the identity: sigma = id with e = +1 is not a quotient
Append(~controls, <"identity map",     DiagonalMatrix(Q, [1,1,1])>);
//  (4) right sigma, wrong y-scale: e^2 f /= f(sigma)*(cX+d)^(2g+2)
Append(~controls, <"wrong y-scale e",  DiagonalMatrix(Q, [-1,3,1])>);
//  (5) a scaling that is not weight-respecting on P(1,g+1,1)
Append(~controls, <"bad weight",       DiagonalMatrix(Q, [-1,1,2])>);

for c in controls do
    name, M := Explode(c);
    okc, _, note := QuotientByInvolution(f14, M);
    error if okc,
        Sprintf("NEGATIVE CONTROL DID NOT FIRE: %o was ACCEPTED -- the helper cannot go red", name);
    nControl +:= 1;
end for;
printf "  %o negative control(s) fired\n", nControl;

// =============================================================== PART 2: 14_1, all three branches
// w_2  : sigma = -x, sum degenerates -> u = x*sigma(x);  y descends       -> genus 0
// w_7  : sigma = -x,                                     y anti-invariant -> genus 1
// w_14 : sigma = id, e = -1, the hyperelliptic involution                 -> P^1
expect14 := [* <2,  HyperellipticCurve(-x^2 + 13*x - 128)>,
               <7,  HyperellipticCurve(x*(-x^2 + 13*x - 128))>,
               <14, HyperellipticCurve(x^2 - x)> *];
ws14 := AssociativeArray();
ws14[2]  := DiagonalMatrix(Q, [-1, 1, 1]);
ws14[7]  := DiagonalMatrix(Q, [-1,-1, 1]);
ws14[14] := DiagonalMatrix(Q, [ 1,-1, 1]);

for e in expect14 do
    m, Cexp := Explode(e);
    okd, Fq, note := QuotientByInvolution(f14, ws14[m]);
    error if not okd, Sprintf("14_1: deriving the quotient by w_%o failed: %o", m, note);
    same, why := QuotientMatches(Fq, Cexp);
    error if not same, Sprintf("14_1: quotient by w_%o disagrees with the published curve (%o)", m, why);
    nDeriv +:= 1;
end for;

// =============================================================== PART 3: 6_7, group generation
// X_0(6,7): y^2 = -3x^4-34x^2-2187.  THREE generators must close up to all seven involutions,
// and every one of the six non-trivial quotients must land on the published curve.
// ⚠ w_m * w_n = w_{mn/gcd(m,n)^2}, not w_{mn}.
f67 := -3*x^4 - 34*x^2 - 2187;
gens67 := AssociativeArray();
gens67[42] := DiagonalMatrix(Q, [1,-1,1]);                   // (x, -y)
gens67[3]  := DiagonalMatrix(Q, [-1,1,1]);                   // (-x, y)
gens67[6]  := Matrix(Q, 3,3, [0,0,1,  0,-27,0,  -27,0,0]);   // (-27/x, -27y/x^2)

all67, consistent, why := ALMatrixGroupFromGenerators(gens67, 2);   // y has weight g+1 = 2
error if not consistent, Sprintf("6_7: generated AL group is inconsistent: %o", why);
error if #Keys(all67) ne 7,
    Sprintf("6_7: three generators closed to %o elements, expected the full group of 7",
            #Keys(all67));
error if ALCompose(3, 6) ne 2 or ALCompose(6, 42) ne 7 or ALCompose(3, 42) ne 14,
    "ALCompose is not computing m*n/gcd(m,n)^2";

expect67 := [* <3,  HyperellipticCurve(-3*x^2 - 34*x - 2187)>,       // u = x^2,      y descends
               <6,  HyperellipticCurve(-3*x^2 - 196)>,               // u = x - 27/x, v = y/x
               <21, HyperellipticCurve(-3*x^2 + 128)>,               // u = x + 27/x, v = y/x
               <7,  HyperellipticCurve((-3*x^2 - 196)*(x^2 + 108))>,
               <2,  HyperellipticCurve((-3*x^2 + 128)*(x^2 - 108))>,
               <42, HyperellipticCurve(x^2 - x)> *];                 // P^1
for e in expect67 do
    m, Cexp := Explode(e);
    error if not IsDefined(all67, m), Sprintf("6_7: w_%o was not generated", m);
    okd, Fq, note := QuotientByInvolution(f67, all67[m]);
    error if not okd, Sprintf("6_7: deriving the quotient by w_%o failed: %o", m, note);
    same, why2 := QuotientMatches(Fq, Cexp);
    error if not same, Sprintf("6_7: quotient by w_%o disagrees with the published curve (%o)", m, why2);
    nDeriv +:= 1;
end for;

// =============================================================== count guards
// ⚠ Raise these when the test grows; NEVER lower one to make a run pass.  The failure mode this
// repo keeps hitting is a check that quietly stops checking.
error if nControl lt 5, Sprintf("only %o negative controls fired, expected 5", nControl);
error if nDeriv lt 9,   Sprintf("only %o quotients were derived and compared, expected 9", nDeriv);
printf "  %o quotient(s) derived and matched against published curves\n", nDeriv;
