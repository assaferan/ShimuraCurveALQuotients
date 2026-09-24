// Kudla-Yang's local Eisenstein factors, written out from the paper rather than from the lattice.
//
// Reference: S. Kudla, T. Yang, "Eisenstein series for SL(2)", Sci. China Math. 53 (2010) 2275-2316.
//
// WHY THIS EXISTS.  SchoferFormula.m computes the local factor as a Whittaker function OF THE LATTICE
// (Wpoly = KY Theorem 4.3, Wpoly2 = Theorem 4.4).  That is correct wherever the local quaternion
// algebra is SPLIT.  At a RAMIFIED prime p | D it is not the right object: the lattice Whittaker can
// VANISHE where the Eisenstein coefficient does not -- concretely, on X0^6(5) at m = 10 and on
// X0^15(2) at m = 10 the code's factor at p = 3 is 0 (there ord_p(m) = 0 and p is inert in
// Q(sqrt m), so its value is 1 + chi(p) = 0), while the Eisenstein coefficient solved from the
// vector-valued oracle is -6 and -4 respectively.  Proposition 5.3 is the statement that separates
// the two cases, and the ANISOTROPIC branch is what this file supplies.
//
// PROPOSITION 5.3 (B a quaternion algebra over Q_p, V its trace-zero part, Q(x) = kappa det x):
//     W_p(s + 1/2, m, mu) = char(Q(mu) + Z_p)(m) * L_p(s + 1, chi_{kappa m})
//                            * b_p(kappa m, s + 1; B) * { zeta_p(2s+2)^{-1}  if B is split,
//                                                         1                  if B is ramified }
// so the ramified case differs from the split one by more than a change of b_p: the local zeta
// factor is ABSENT.  Since SchoferFormula.m divides by en/ed = prod L_p^{-1} * prod zeta_p, that
// missing factor is over-corrected at p | D, on top of b_p itself being wrong there.
//
// The b_p are KY (2.16) and (2.17).  With 4*kappa*m = d*c^2 (eq (2.15), d the fundamental
// discriminant of Q(sqrt(kappa m))), k = ord_p(c), X = p^{-s}, and v_p = 1, 0, -1 according as
// Q(sqrt(kappa m))/Q is split, ramified or inert at p:
//   p does not divide D:
//     b_p = [ 1 - v X + p^k v X^{1+2k} - p^{k+1} X^{2k+2} ] / (1 - p X^2)
//   p divides D:
//     b_p = [ (1 - v X)(1 - p^2 X^2) - v p^{k+1} X^{2k+1} + p^{k+2} X^{2k+2}
//                                    + v p^{k+1} X^{2k+3} - p^{2k+2} X^{2k+4} ] / (1 - p X^2)
//
// CONVENTIONS, and they are easy to get wrong (see tests/KudlaYangLocal.m):
//   * the character and the conductor are those of Q(sqrt(kappa m)).  The pipeline passes negQ = -Q
//     to its Whittaker routine, which makes the relevant field Q(sqrt m) -- NOT the Q(sqrt(-m)) used
//     a few lines away for the class number h/w.  These intrinsics take kappa*m directly, so the
//     caller must pass the value whose square root generates the field it means.
//   * k is the exponent of the CONDUCTOR, not ord_p(m): m = 25 at p = 5 has d = -4, c = 5, k = 1.
//   * the weight-3/2 special value is s = 0, so W_p sits at X = 1 while b_p, whose argument is s + 1,
//     sits at X = 1/p.
//
// WHAT SWAPPING IN THE RAMIFIED BRANCH DOES *NOT* FIX -- record it, it is informative.
// Replacing the lattice Whittaker at p | D by Prop 5.3's anisotropic value does NOT reproduce the
// Eisenstein coefficients solved from the vector-valued oracle (memory b-eisenstein-coefficients-
// solved).  With kappa*m = +r the ramified factor still vanishes on X0^6(5) at r = 10 (p = 3 is
// SPLIT in Q(sqrt 10) with conductor exponent 0, and (2.17) vanishes there), while the true
// coefficient is -6.  With kappa*m = -r that case is fixed but X0^15(2) at r = 2 breaks instead.
// No single sign convention works, and no variant of the assembly is proportional to the true
// coefficients.  So the residual defect is NOT merely "the wrong local factor at p | D": a product of
// local densities of ONE quadratic space cannot produce them.
// The natural explanation, and the next thing to test: for signature (1,2) the obstruction Eisenstein
// series is INCOHERENT -- attached to a collection differing from a coherent one at a single place --
// so it vanishes at the centre and its coefficients carry a DERIVATIVE W'_p at the distinguished
// place rather than a plain product of values.  That also fits the standing observation that the
// pipeline computes the Whittaker polynomial and then discards its s-derivative.

// Fundamental discriminant of Q(sqrt(a)) for a rational a, and the c of 4a = d c^2.
function disc_and_conductor(a)
    sq := SquarefreeFactorization(Numerator(a) * Denominator(a));
    d := (sq eq 1) select 1 else Discriminant(Integers(QuadraticField(sq)));
    is_sq, c := IsSquare(Rationals()!(4*a) / d);
    error if not is_sq, Sprintf("4a/d is not a square for a = %o (d = %o)", a, d);
    return d, c;
end function;

intrinsic KYSplittingType(kappam::FldRatElt, p::RngIntElt) -> RngIntElt
{v_p in 1, 0, -1 according as Q(sqrt(kappa m))/Q is split, ramified or inert at p (KY Section 2).}
    require IsPrime(p) : "p must be prime.";
    d := disc_and_conductor(kappam);
    if d eq 1 then return 1; end if;                  // Q(sqrt 1) = Q, split everywhere
    return KroneckerSymbol(d, p);
end intrinsic;

intrinsic KYbp(kappam::FldRatElt, p::RngIntElt, ramified::BoolElt, X::FldRatElt) -> FldRatElt
{Kudla-Yang's local density b_p(kappa m, s; D) at X = p^(-s), from (2.16) when p does not divide D
 and (2.17) when it does.  `ramified` says whether the quaternion algebra is ramified at p, i.e.
 whether p divides D.}
    require IsPrime(p) : "p must be prime.";
    require kappam ne 0 : "kappa*m must be nonzero.";
    d, c := disc_and_conductor(kappam);
    k := Valuation(c, p);
    require k ge 0 : "KY (2.16)/(2.17) assume ord_p(c) >= 0.";
    v := (d eq 1) select 1 else KroneckerSymbol(d, p);
    P := Rationals()!p;
    den := 1 - P*X^2;
    require den ne 0 : "b_p has a pole at this X.";
    if not ramified then
        num := 1 - v*X + P^k*v*X^(1+2*k) - P^(k+1)*X^(2*k+2);
    else
        num := (1 - v*X)*(1 - P^2*X^2)
               - v*P^(k+1)*X^(2*k+1) + P^(k+2)*X^(2*k+2)
               + v*P^(k+1)*X^(2*k+3) - P^(2*k+2)*X^(2*k+4);
    end if;
    return num/den;
end intrinsic;

intrinsic KYWhittaker53(kappam::FldRatElt, p::RngIntElt, ramified::BoolElt : s := 0) -> FldRatElt
{The local Whittaker value W_p(s + 1/2, m, mu) of Kudla-Yang Proposition 5.3, for mu with
 m in Q(mu) + Z_p (so the characteristic function is 1) -- in particular for mu = 0 and integral m.
 The default s = 0 is the weight-3/2 special value the CM-value formula needs.}
    require IsPrime(p) : "p must be prime.";
    P := Rationals()!p;
    Xb := P^(-(s+1));                                   // b_p's argument is s + 1
    bp := KYbp(kappam, p, ramified, Xb);
    d := disc_and_conductor(kappam);
    chi := (d eq 1) select 1 else KroneckerSymbol(d, p);
    Lp := 1/(1 - chi*P^(-(s+1)));                       // L_p(s+1, chi_{kappa m})
    if ramified then
        return Lp * bp;
    end if;
    return Lp * bp * (1 - P^(-(2*s+2)));                // times zeta_p(2s+2)^{-1}
end intrinsic;

intrinsic KYWhittaker54(kappam::FldRatElt, p::RngIntElt, e::RngIntElt) -> FldRatElt
{Kudla-Yang Proposition 5.4 at X = 1: the local Whittaker at a prime where the quaternion algebra is
 SPLIT but the lattice has Eichler level p^e (so p | N).  Written for e = 1, the squarefree-level case
 this project needs.  With 4*kappa*m = d*c^2 and k = ord_p(c),
     p | d      : 1 + (1-1/p) sum over 1<=l<=k of p^(1-l), minus p^(-k-1)
     p not | d  : 1 + (1-1/p) sum over 1<=l<=k of p^(1-l), plus chi_d(p) p^(-k)
 which for k = 0 collapses to 1 - 1/p and 1 + chi_d(p).  Cross-checked against the pipeline's own
 factor by tests/KudlaYangLocal.m.}
    require IsPrime(p) : "p must be prime.";
    require e eq 1 : "Only the Eichler level e = 1 (squarefree N) is implemented.";
    d, c := disc_and_conductor(kappam);
    k := Valuation(c, p);
    require k ge 0 : "Assumes ord_p(c) >= 0.";
    chi := (d eq 1) select 1 else KroneckerSymbol(d, p);
    P := Rationals()!p;
    total := Rationals()!1;
    for l in [1..k] do total +:= (1 - 1/P) * P^(1-l); end for;
    if (Integers()!d) mod p eq 0 then
        return total - P^(-k-1);
    end if;
    return total + chi * P^(-k);
end intrinsic;
