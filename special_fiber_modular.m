// Special-fiber (reduction mod p) test for non-hyperellipticity of AL quotients of
// modular curves, generalizing Furumoto-Hasegawa, "Hyperelliptic Quotients of Modular
// Curves", Section 5.  Currently for D = 1 (classical modular curves) with a genus-0
// special-fiber component.
//
// The conclusion is GEOMETRIC non-hyperellipticity (not just over Q).  For genus >= 2 the
// hyperelliptic involution is unique, hence Galois-stable, hence defined over Q regardless
// of whether the curve is hyperelliptic over Q or only over Qbar; it therefore reduces mod p
// to an F_p-rational automorphism inducing an F_p-rational, Frobenius-compatible Mobius
// involution on the genus-0 component.  Ruling out the existence of such an involution (below)
// rules out the reduction of any geometric hyperelliptic involution.  Cross-checked for D = 1
// against the canonical-model (Petri) oracle in verify_d1_hyperelliptic.m, which likewise
// decides geometric hyperellipticity; consistent with the geometric IsHyp flag used pipeline-
// wide (see ShimuraQuotients.m, "Committed to Geometrically Hyperelliptic curves").
//
// Method.  Fix a prime p with p | N.  X_0(N)/W reduces mod p to a curve whose components
// are copies of X_0(N/p)/W''  where  W'' = { w / p^v_p(w) : w in W }  is the p-coprime
// image of W (the component group).  Two cases:
//
//   Case 1 (no w in W is divisible by p):  the special fiber is two components, each
//     isomorphic to X_0(N/p)/W'', crossing at the supersingular points.  If X_0(N)/W is
//     hyperelliptic with involution u, then g = w_p * u fixes each component and acts on
//     it as Frobenius on the supersingular points: it fixes the F_p-rational ones and
//     conjugates the others.  ALL supersingular points are constrained.
//
//   Case 2 (some w in W is divisible by p):  the special fiber is a single component
//     X_0(N/p)/W'', self-intersecting at the supersingular points.  Here the hyperelliptic
//     involution u itself acts as Frobenius on the supersingular points that are NOT
//     F_p-rational (the F_p-rational ones are free).  Only the non-F_p ss are constrained.
//
// In both cases the action on the genus-0 component is a Mobius transformation g, pinned
// down by its values on >= 3 of the constrained points.  If the unique such g fails to send
// every constrained point to its Frobenius image, no compatible involution exists and the
// curve is NOT geometrically hyperelliptic.  If g comes out as the identity (Frobenius acts
// trivially on the constrained points) the test is inconclusive; in case 1 g = id is the
// admissible possibility u = w_p.  A non-F_p-rational g also yields "not geometrically
// hyperelliptic" (the geometric hyperelliptic involution would reduce to an F_p-rational g).
//
// Validation (data/curves_after_UpdateCurves8.dat, D = 1, genus >= 3): reproduces FH on
// X_0(42)/<w_2> (case 1) and X_0(42)/<w_7> (case 2); 0 contradictions against the 59
// known-hyperelliptic curves; proves 14 of the 43 then-undetermined curves non-hyperelliptic,
// every one independently confirmed by verify_d1_hyperelliptic.m (canonical-model oracle).

// ---- P^1 helpers ----
// Mobius matrix sending the standard frame to the ordered triple t.
function aut01oo(t)
  p0:=t[1]; p1:=t[2]; poo:=t[3];
  c0:=p1[1]*p0[2]-p0[1]*p1[2]; d0:=poo[1]*p1[2]-p1[1]*poo[2];
  return Matrix([[c0*poo[1], p0[1]*d0],[c0*poo[2], d0*p0[2]]]);
end function;
// Mobius matrix sending ordered triple o to ordered triple d.
function mobm(o,d) return aut01oo(d)*aut01oo(o)^(-1); end function;
// apply Mobius matrix g to a projective point pt = [x:y].
app := func<g,pt | [g[1,1]*pt[1]+g[1,2]*pt[2], g[2,1]*pt[1]+g[2,2]*pt[2]]>;
// projective equality of two points on P^1.
eqP := func<a,b | a[1]*b[2]-a[2]*b[1] eq 0>;

// supersingular points of X_0(M) over F_{p^2} (the special-fiber crossing points).
function ss_points(M, p, F)
  X := SmallModularCurve(M); j := jInvariant(X,M);
  Xp := ChangeRing(X,F); RR := CoordinateRing(Ambient(Xp));
  jn := RR!Numerator(j); jd := RR!Denominator(j);
  ss := [];
  for j0 in [r[1] : r in Roots(ChangeRing(SupersingularPolynomial(p),F))] do
    for pt in Points(Scheme(Xp, jn-j0*jd)) do e := Eltseq(pt);
      if not exists{q:q in ss|eqP(e,q)} then Append(~ss, e); end if; end for; end for;
  return ss;
end function;

// push ss points through the quotient by a single AL involution w_q on X_0(M).
function push_single(ss, M, q, F)
  X := SmallModularCurve(M); R := CoordinateRing(Ambient(X));
  dp := DefiningPolynomials(AtkinLehnerInvolution(X, M, q));
  a:=F!MonomialCoefficient(dp[1],R.1); b:=F!MonomialCoefficient(dp[1],R.2);
  c:=F!MonomialCoefficient(dp[2],R.1); d:=F!MonomialCoefficient(dp[2],R.2);
  qss := [];
  for P in ss do Q := [c*P[1]^2+(a+d)*P[1]*P[2]+b*P[2]^2, c*P[1]*P[2]+d*P[2]^2];
    if Q ne [F!0,F!0] and not exists{T:T in qss|eqP(Q,T)} then Append(~qss,Q); end if; end for;
  return qss;
end function;

// push ss points through the quotient by the whole component group W'' (order 1, 2 or 4) on
// X_0(M).  For |W''| = 2 this is the single-involution quotient; for the Klein-four case |W''| = 4
// we quotient by one generator, then by the second generator's INDUCED involution on the first
// quotient (the AL involutions commute, so w_q2 descends; its matrix is fitted from 3 image pairs).
function push_group(ss, M, Wpp, F)
  if #Wpp eq 1 then return ss; end if;
  X := SmallModularCurve(M); R := CoordinateRing(Ambient(X));
  almat := function(q)
    dp := DefiningPolynomials(AtkinLehnerInvolution(X, M, q));
    return Matrix(F,2,2,[MonomialCoefficient(dp[1],R.1),MonomialCoefficient(dp[1],R.2),
                         MonomialCoefficient(dp[2],R.1),MonomialCoefficient(dp[2],R.2)]);
  end function;
  // quotient X_0(M) -> X_0(M)/<involution m> on P^1 coordinates:
  qmap := func< m,P | [m[2,1]*P[1]^2+(m[1,1]+m[2,2])*P[1]*P[2]+m[1,2]*P[2]^2, m[2,1]*P[1]*P[2]+m[2,2]*P[2]^2] >;
  gens := Sort([q : q in Setseq(Wpp) | q ne 1]);
  m1 := almat(gens[1]);
  cur := [];
  for P in ss do Q := qmap(m1,P);
    if Q ne [F!0,F!0] and not exists{T:T in cur|eqP(Q,T)} then Append(~cur,Q); end if; end for;
  if #Wpp eq 2 then return cur; end if;
  // |W''| = 4 : induced second involution h on X_0(M)/<gens[1]>, fitted from 3 image pairs
  m2 := almat(gens[2]);
  pairs := [];
  for P in ss do sv := qmap(m1,P); dv := qmap(m1, app(m2,P));
    if sv ne [F!0,F!0] and dv ne [F!0,F!0] and not exists{pr:pr in pairs|eqP(pr[1],sv)} then
      Append(~pairs, <sv,dv>); end if; end for;
  assert #pairs ge 3;
  h := mobm([pairs[1][1],pairs[2][1],pairs[3][1]], [pairs[1][2],pairs[2][2],pairs[3][2]]);
  cur2 := [];
  for Q in cur do QQ := qmap(h,Q);
    if QQ ne [F!0,F!0] and not exists{T:T in cur2|eqP(QQ,T)} then Append(~cur2,QQ); end if; end for;
  return cur2;
end function;

// Does a compatible (order-2, rational, = Frobenius-on-ss) Mobius involution exist?
// Returns "yes" (compatible -> inconclusive), "no" (none -> NOT hyperelliptic), or
// "underdet" (fewer than 3 constrained points, or g = identity -> inconclusive).
function has_compat(qss, p, F, case2)
  frob := func<pt | [Frobenius(pt[1]),Frobenius(pt[2])]>;
  pts := case2 select [Q : Q in qss | not eqP(frob(Q),Q)] else qss;   // case 2: non-F_p ss only
  if #pts lt 3 then return "underdet"; end if;
  g := mobm([pts[1],pts[2],pts[3]], [frob(pts[1]),frob(pts[2]),frob(pts[3])]);
  if not &and[eqP(app(g,Q),frob(Q)) : Q in pts] then return "no"; end if;
  if g[1,2] eq 0 and g[2,1] eq 0 and g[1,1] eq g[2,2] then return "underdet"; end if;
  s := exists(c){x : x in [g[1,1],g[1,2],g[2,1],g[2,2]] | x ne 0} select c else F!1;
  if not &and[x in GF(p) : x in [g[1,1]/s,g[1,2]/s,g[2,1]/s,g[2,2]/s]] then return "no"; end if;
  return "yes";
end function;

intrinsic SpecialFiberNotHyperelliptic(N::RngIntElt, W::SetEnum) -> BoolElt, MonStgElt
{Special-fiber reduction-mod-p test for X_0(N)/W (D = 1).  Returns true with a witness
string if the curve is proven NOT geometrically hyperelliptic (a fortiori not hyperelliptic
over Q), otherwise false.  Reaches any prime p | N (including composite N) for which X_0(N/p) is
genus 0 and the special-fiber component X_0(N/p)/W'' is genus 0; the component group W'' may be
trivial, an involution, or a Klein-four group (order up to 4).  Returns false otherwise.  (The
case where X_0(N/p) has positive genus but X_0(N/p)/W'' is genus 0 needs supersingular points on
the quotient and is a TODO.)}
    for p in PrimeDivisors(N) do
        M := N div p;
        if GenusShimuraCurveQuotient(1, M, {Integers()|1}) ne 0 then continue; end if;   // need X_0(M) = P^1
        Wpp := {w div p^Valuation(w,p) : w in W};            // p-coprime image of W = component group
        if GenusShimuraCurveQuotient(1, M, Wpp) ne 0 then continue; end if;
        if #Wpp gt 4 then continue; end if;                  // larger group quotients: TODO
        case2 := exists{w : w in W | w mod p eq 0};
        F := GF(p^2);
        ss := ss_points(M, p, F);
        qss := push_group(ss, M, Wpp, F);
        if has_compat(qss, p, F, case2) eq "no" then
            return true, Sprintf("SpecialFiber p=%o case=%o component=X_0(%o)/%o", p, case2 select 2 else 1, M, Wpp);
        end if;
    end for;
    return false, _;
end intrinsic;

intrinsic FilterBySpecialFiber(~curves::SeqEnum)
{Mark D = 1, 6, 10 and 22 AL quotients proven (geometrically) non-hyperelliptic by the
special-fiber reduction-mod-p test (Furumoto-Hasegawa Section 5, generalized).  D = 1 uses
SpecialFiberNotHyperelliptic (classical modular curves), D = 6 uses SpecialFiberNotHyperellipticD6
(hypergeometric), D = 10 uses SpecialFiberNotHyperellipticD10 (Heun), D = 22 uses
SpecialFiberNotHyperellipticD22 (general CM-value, coverage-limited); these handle the base level
M = 1.  For any D, the generic SpecialFiberNotHyperellipticCM then handles level-M > 1 components
X_0(D,M)/W'' (Phase B) for bases registered in CM_BASE_DATA.  Only curves with a genus-0 special-fiber
component are reached, everything else is left untouched.}
    for i->X in curves do
        if assigned X`IsSubhyp then continue; end if;
        if X`g lt 3 then continue; end if;
        if X`D eq 1 then
            ok, witness := SpecialFiberNotHyperelliptic(X`N, X`W);
        elif X`D eq 6 then
            ok, witness := SpecialFiberNotHyperellipticD6(X`N, X`W);
        elif X`D eq 10 then
            ok, witness := SpecialFiberNotHyperellipticD10(X`N, X`W);
        elif X`D eq 22 then
            ok, witness := SpecialFiberNotHyperellipticD22(X`N, X`W);
        else
            ok := false;
        end if;
        if not ok then                           // generic level-M > 1 CM-value test (Phase B)
            ok, witness := SpecialFiberNotHyperellipticCM(X`D, X`N, X`W);
        end if;
        if ok then
            curves[i]`IsSubhyp := false;
            curves[i]`IsHyp := false;
            curves[i]`TestInWhichProved := witness;
        end if;
        if (i mod 200 eq 0) then
            vprintf ShimuraQuotients, 1: "i = %o/%o\n", i, #curves;
        end if;
    end for;
end intrinsic;

// ---------------------------------------------------------------------------
// D = 6 (quaternionic) special-fiber test.  The supersingular points of the
// genus-0 Shimura curve X_0(6,1) at a prime p (p not dividing 6) are computed
// via Elkies, "Shimura Curve Computations" Section 3.3: the [p/24] roots of the
// (2,4,6) hypergeometric polynomial give the generic supersingular points on the
// star X_0^*(6,1), and the three elliptic vertices tau = 0,1,oo (CM discriminants
// -24, -4, -3) are supersingular iff p is inert in the corresponding quadratic
// order.  On X_0(6,1) (the degree-4 Klein-four cover of the star) each generic
// star point lifts to its full 4-point orbit via (x:y:z) -> (x^2:y^2:z^2), i.e.
// z^2 = -tau x^2, and each supersingular vertex to its 2-point orbit.  The total
// is asserted to equal genus(X_0(6,p)) + 1 (the number of nodes of the two-
// component special fiber), which is a complete arithmetic self-check.
//
// The Mobius test runs on the W'' = {1} component (X_0(6,1), the conic x^2+3y^2+z^2),
// on the full-Atkin-Lehner star component (the tau-line), and on the three intermediate
// quotients X_0(6,1)/<w_q>, q in {2,3,6}.  Each w_q acts on the conic as a diagonal sign
// flip (w_2 flips y, w_3 flips x, w_6 flips x and y == flips z projectively), so its two
// linear invariant coordinates give the quotient map conic -> P^1 directly: w_2 -> (x:z),
// w_3 -> (y:z), w_6 -> (x:y).  The pushed supersingular-point count is checked against
// genus(X_0(6,p)/<w_q>) + 1 on every call (the same self-check as the conic case).
// As in the D=1 case the conclusion is GEOMETRIC non-hyperellipticity: the (unique, hence
// Q-rational) hyperelliptic involution would reduce to an F_p-rational Frobenius-compatible
// Mobius involution, and the test rules such out.  For D=6 this is intrinsic -- the genus-0
// components are pointless conics over Q (x^2+3y^2+z^2 has no Q-point), yet the test works on
// their reduction mod p where every conic acquires a rational point, so it cannot be sensing
// a mere Q-rationality obstruction.
//
// Validated on data/curves_after_UpdateCurves8.dat (D=6, genus >= 3): proves X_0(6,23)/<w_23>
// non-hyperelliptic (conic case); the intermediate-quotient branch proves 112 curves, with 0
// contradictions against the 59 reached geometrically-hyperelliptic curves (zero false
// positives), and newly settles X_0(6,229)/<w_6,w_229>.

d6_inert := func<p,d | not IsSplit(p, MaximalOrder(QuadraticField(d)))>;

function d6_hypg(p, F)
  if   p mod 24 in [19,23] then A:=19/24; B:=23/24; C:=3/2;
  elif p mod 24 in [1,5]   then A:=1/24;  B:=5/24;  C:=1/2;
  elif p mod 24 in [7,11]  then A:=7/24;  B:=11/24; C:=1/2;
  else                          A:=13/24; B:=17/24; C:=3/2; end if;
  R := PolynomialRing(F);
  return &+[ &*[F | (A+l)*(B+l)/(C+l) : l in [0..j-1]]/Factorial(j) * R.1^j : j in [0..(p div 24)]];
end function;

// supersingular points of X_0(6,1) on the conic x^2+3y^2+z^2 over F (= F_{p^2}).
function d6_ss_conic(p, F)
  eqP3 := func<a,b | a[1]*b[2]-a[2]*b[1] eq 0 and a[1]*b[3]-a[3]*b[1] eq 0 and a[2]*b[3]-a[3]*b[2] eq 0>;
  pts := []; add := procedure(~pts,P) if not exists{q:q in pts|eqP3(P,q)} then Append(~pts,P); end if; end procedure;
  for r in Roots(d6_hypg(p,F)) do tau0:=r[1]; rz:=Sqrt(-tau0); ry:=Sqrt(-(1-tau0)/3);
    for sy in [ry,-ry] do for sz in [rz,-rz] do add(~pts,[F!1,sy,sz]); end for; end for; end for;
  if d6_inert(p,-3)  then s:=Sqrt(F!-3); add(~pts,[F!0,F!1,s]); add(~pts,[F!0,F!1,-s]); end if;
  if d6_inert(p,-4)  then s:=Sqrt(F!-1); add(~pts,[s,F!0,F!1]); add(~pts,[-s,F!0,F!1]); end if;
  if d6_inert(p,-24) then s:=Sqrt(F!-3); add(~pts,[s,F!1,F!0]); add(~pts,[-s,F!1,F!0]); end if;
  return pts;
end function;

// push conic SS points to an intermediate quotient X_0(6,1)/<w_q>, q in {2,3,6}.  w_q is a
// diagonal sign flip, so the quotient map to P^1 is the pair of linear invariant coordinates:
// w_2 (flip y) -> (x:z), w_3 (flip x) -> (y:z), w_6 (flip x,y) -> (x:y).
function d6_proj_quot(ss, q, F)
  idx := case< q | 2: [1,3], 3: [2,3], 6: [1,2], default: [Integers()|] >;
  qss := [];
  for P in ss do Q := [P[idx[1]], P[idx[2]]];
    if Q ne [F!0,F!0] and not exists{T:T in qss|eqP(Q,T)} then Append(~qss,Q); end if; end for;
  return qss;
end function;

// supersingular points of the star X_0^*(6,1) in the hypergeometric tau-coordinate (P^1).
function d6_ss_star(p, F)
  pts := [[r[1],F!1] : r in Roots(d6_hypg(p,F))];
  if d6_inert(p,-24) then Append(~pts,[F!0,F!1]); end if;
  if d6_inert(p,-4)  then Append(~pts,[F!1,F!1]); end if;
  if d6_inert(p,-3)  then Append(~pts,[F!1,F!0]); end if;
  return pts;
end function;

// project the conic SS points to P^1 over an F_p-rational point (the conic has no Q-point).
function d6_proj_P1(p, ssF)
  Fp := GF(p); F := GF(p^2); P2<x,y,z> := ProjectiveSpace(Rationals(),2); Cm := Curve(P2, x^2+3*y^2+z^2);
  C := ChangeRing(Cm, Fp); Cc,CtoCc := Conic(C); _,P0 := HasRationalPoint(Cc); CP1,pi := Projection(Cc,P0);
  C2 := ChangeRing(Cm,F); C2c := ChangeRing(Cc,F); C2P1 := ChangeRing(CP1,F);
  am := AlgebraMap(CtoCc); m1 := map<C2->C2c | [am(Domain(am).i):i in [1..3]]>;
  am := AlgebraMap(pi);     m2 := map<C2c->C2P1 | [am(Domain(am).i):i in [1..2]]>;
  fmap := m1*m2;
  return [Eltseq(fmap(C2!Eltseq(P))) : P in ssF];
end function;

intrinsic SpecialFiberNotHyperellipticD6(N::RngIntElt, W::SetEnum) -> BoolElt, MonStgElt
{Special-fiber reduction-mod-p test for the discriminant-6 Shimura curve quotient
X_0(6,N)/W.  Returns true with a witness if proven NOT geometrically hyperelliptic
(a fortiori not hyperelliptic over Q), otherwise false.  Implemented for N = p prime
with p-coprime component group W'' equal to the
trivial group, an intermediate single-involution quotient by w_q (q = 2, 3 or 6), or
the full Atkin-Lehner group; other cases return false.  Except for the full-star case
the pushed supersingular-point count is checked against genus of X_0(6,p)/w_q plus 1
on every call.}
    for p in PrimeDivisors(N) do
        if 6 mod p eq 0 or N div p ne 1 then continue; end if;
        Wpp := {w div p^Valuation(w,p) : w in W};
        case2 := exists{w : w in W | w mod p eq 0};
        F := GF(p^2);
        if Wpp eq {Integers()|1} then
            ss := d6_ss_conic(p, F);
            assert #ss eq GenusShimuraCurveQuotient(6, p, {Integers()|1}) + 1;   // arithmetic self-check
            qss := d6_proj_P1(p, ss);
        elif Wpp in {{Integers()|1,2},{Integers()|1,3},{Integers()|1,6}} then
            ss := d6_ss_conic(p, F);
            qss := d6_proj_quot(ss, Rep(Wpp diff {Integers()|1}), F);
            assert #qss eq GenusShimuraCurveQuotient(6, p, Wpp) + 1;             // arithmetic self-check
        elif Wpp eq {1,2,3,6} then
            qss := d6_ss_star(p, F);
        else
            continue;                            // remaining cases (composite N etc.): TODO
        end if;
        if has_compat(qss, p, F, case2) eq "no" then
            return true, Sprintf("SpecialFiberD6 p=%o case=%o component=X_0(6,1)/%o", p, case2 select 2 else 1, Wpp);
        end if;
    end for;
    return false, _;
end intrinsic;

// ---------------------------------------------------------------------------
// Quaternionic (Brandt-module) supersingular layer, general discriminant D.
//
// For a quaternionic discriminant D (squarefree, even number of prime factors) and a prime
// p not dividing D, the supersingular points of the Shimura curve X_0(D,1) in characteristic
// p are the left-ideal classes of a maximal order in the DEFINITE quaternion algebra of
// discriminant D*p (ramified at the primes of D, at p, and at infinity).  The Brandt module
// B(D*p) realises this:
//
//   * dim B(D*p) = #supersingular points = genus(X_0(D,p)) + 1  (the node count of the two-
//     component special fiber);
//   * the geometric Frobenius at p acts on the supersingular points as the Atkin-Lehner
//     involution w_p of B(D*p)  (Deuring/Eichler);
//   * each Atkin-Lehner involution w_q of X_0(D,1) (q | D a Hall divisor) acts as the w_q of
//     B(D*p).
//
// Magma's AtkinLehnerOperator on the Brandt module returns, at each ramified prime power, a
// signed permutation matrix (uniform sign -1 here); the underlying point permutation is read
// off from the support, so the sign is irrelevant for the permutation itself (it only flips
// the trace).  Composite involutions (w_6, w_10, ...) are products of the prime ones.
//
// Validated against the D=6 geometric layer (d6_ss_conic + sign-flip AL involutions): for all
// primes 7..59 the Brandt fixed-point counts of w_2, w_3, w_6 AND of Frobenius agree exactly
// with the geometric counts.  This layer supplies the supersingular SET with its Frobenius and
// Atkin-Lehner action for discriminants beyond D=6 (e.g. D=10, whose star is the (2,2,2,3)
// quadrilateral orbifold and so has no Elkies hypergeometric model); the remaining ingredient
// for a full non-hyperellipticity test -- the points' coordinates on a genus-0 model -- is a
// separate (harder) step.

// permutation in Sym(n) carried by a permutation matrix Mx (one nonzero entry per row).
function ss_perm_of_matrix(Mx, n)
  return Sym(n) ! [ [j : j in [1..n] | Mx[i][j] ne 0][1] : i in [1..n] ];
end function;

intrinsic SupersingularALData(D::RngIntElt, p::RngIntElt) -> RngIntElt, GrpPermElt, Assoc
{Supersingular points of the discriminant-D Shimura curve X_0(D,1) in characteristic p,
computed via the Brandt module of the definite quaternion algebra of discriminant D*p
(p prime, p not dividing D; D squarefree with an even number of prime factors).  Returns:
(1) the number n of supersingular points; (2) the geometric Frobenius at p as a permutation
in Sym(n); (3) an associative array sending each Hall divisor q | D, q > 1, to the
Atkin-Lehner involution w_q as a permutation in Sym(n).  The count is checked against
genus(X_0(D,p)) + 1, and the involution/commutation relations are asserted, on every call.}
    require IsPrime(p): "p must be prime";
    require D gt 1 and IsSquarefree(D) and IsEven(#PrimeDivisors(D)):
        "D must be a quaternionic discriminant (squarefree, even number of primes)";
    require D mod p ne 0: "p must not divide D";

    M := BrandtModule(D*p);
    n := Dimension(M);
    assert n eq GenusShimuraCurveQuotient(D, p, {Integers()|1}) + 1;        // arithmetic self-check

    // prime building blocks: Frobenius = w_p, and w_q for each prime q | D
    primeperm := AssociativeArray();
    for q in PrimeDivisors(D) cat [p] do
        primeperm[q] := ss_perm_of_matrix(Matrix(AtkinLehnerOperator(M, q)), n);
        assert primeperm[q]^2 eq Id(Sym(n));                                // each w_q is an involution
    end for;
    frob := primeperm[p];

    // Atkin-Lehner involutions of X_0(D,1): w_q for every Hall divisor q | D, q > 1.
    al := AssociativeArray();
    for q in [d : d in Divisors(D) | d gt 1] do
        al[q] := &*[Sym(n) | primeperm[r] : r in PrimeDivisors(q)];          // composite = product of primes
        assert al[q]^2 eq Id(Sym(n));
    end for;

    // the curve is defined over Q, so its Atkin-Lehner involutions commute with Frobenius,
    // and the Atkin-Lehner group is abelian.
    for q in Keys(al) do
        assert al[q]*frob eq frob*al[q];
        for r in Keys(al) do assert al[q]*al[r] eq al[r]*al[q]; end for;
    end for;

    return n, frob, al;
end intrinsic;
