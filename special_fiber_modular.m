// Special-fiber (reduction mod p) test for non-hyperellipticity of AL quotients of
// modular curves, generalizing Furumoto-Hasegawa, "Hyperelliptic Quotients of Modular
// Curves", Section 5.  Currently for D = 1 (classical modular curves) with a genus-0
// special-fiber component.  NOT yet wired into the pipeline.
//
// Method.  Fix a prime p with p | N.  X_0(N)/W reduces mod p to a curve whose components
// are copies of X_0(N/p)/W''  where  W'' = { w / p^{v_p(w)} : w in W }  is the p-coprime
// image of W (the component group).  Two cases:
//
//   Case 1 (no w in W is divisible by p):  the special fiber is two components, each
//     isomorphic to X_0(N/p)/W'', crossing at the supersingular points.  If X_0(N)/W is
//     hyperelliptic with involution u, then g = w_p * u fixes each component and acts on
//     it as Frobenius on the supersingular points: it fixes the F_p-rational ones and
//     conjugates the others.  We constrain ALL supersingular points.
//
//   Case 2 (some w in W is divisible by p):  the special fiber is a single component
//     X_0(N/p)/W'', self-intersecting at the supersingular points.  Here the hyperelliptic
//     involution u itself acts as Frobenius on the supersingular points that are NOT
//     F_p-rational (the F_p-rational ones are left free).  We constrain only the non-F_p ss.
//
// In both cases the action on the genus-0 component is a Mobius transformation g, pinned
// down by its values on >= 3 of the constrained points.  If the unique such g fails to
// send every constrained point to its Frobenius image, no compatible involution exists and
// the curve is NOT hyperelliptic.  If g comes out as the identity (Frobenius acts trivially
// on the constrained points) the test is inconclusive; in case 1 g = id corresponds to the
// admissible possibility u = w_p.  A non-rational g also yields "not hyperelliptic" (the
// canonical involution is defined over the base field).
//
// Validation (data/curves_after_UpdateCurves8.dat, D = 1, genus >= 3):
//   * reproduces FH on X_0(42)/<w_2> (case 1) and X_0(42)/<w_7> (case 2);
//   * 0 contradictions against the 59 known-hyperelliptic curves;
//   * proves 14 of the 43 then-undetermined curves non-hyperelliptic, every one of which
//     was independently confirmed by verify_d1_hyperelliptic.m (canonical-model oracle).

SetQuitOnError(true);
AttachSpec("ShimuraQuotients.spec");

// ---- P^1 helpers ----
// Mobius matrix sending the ordered triple (0,1,oo)-preimages; aut01oo(t) maps 0,1,oo to t.
aut01oo := function(t)
  p0:=t[1]; p1:=t[2]; poo:=t[3];
  c0:=p1[1]*p0[2]-p0[1]*p1[2]; d0:=poo[1]*p1[2]-p1[1]*poo[2];
  return Matrix([[c0*poo[1], p0[1]*d0],[c0*poo[2], d0*p0[2]]]);
end function;
// Mobius matrix sending ordered triple o to ordered triple d.
mobm := function(o,d) return aut01oo(d)*aut01oo(o)^(-1); end function;
// apply Mobius matrix g to a projective point pt = [x:y].
app := func<g,pt | [g[1,1]*pt[1]+g[1,2]*pt[2], g[2,1]*pt[1]+g[2,2]*pt[2]]>;
// projective equality of two points on P^1.
eqP := func<a,b | a[1]*b[2]-a[2]*b[1] eq 0>;

// supersingular points of X_0(M) over F_{p^2} (the special-fiber crossing points).
ss_points := function(M, p, F)
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
push_single := function(ss, M, q, F)
  X := SmallModularCurve(M); R := CoordinateRing(Ambient(X));
  dp := DefiningPolynomials(AtkinLehnerInvolution(X, M, q));
  a:=F!MonomialCoefficient(dp[1],R.1); b:=F!MonomialCoefficient(dp[1],R.2);
  c:=F!MonomialCoefficient(dp[2],R.1); d:=F!MonomialCoefficient(dp[2],R.2);
  qss := [];
  for P in ss do Q := [c*P[1]^2+(a+d)*P[1]*P[2]+b*P[2]^2, c*P[1]*P[2]+d*P[2]^2];
    if Q ne [F!0,F!0] and not exists{T:T in qss|eqP(Q,T)} then Append(~qss,Q); end if; end for;
  return qss;
end function;

// Does a compatible (order-2, rational, = Frobenius-on-ss) Mobius involution exist?
// Returns "yes" (compatible -> inconclusive), "no" (none -> NOT hyperelliptic), or
// "underdet" (fewer than 3 constrained points, or g = identity -> inconclusive).
has_compat := function(qss, p, F, case2)
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

GENUS0 := [M : M in [1..50] | GenusShimuraCurveQuotient(1,M,{1}) eq 0];

// Special-fiber test for X_0(N)/W (D = 1).  Returns true with a witness tuple
// <p, case, M, W''> if proven NOT hyperelliptic; otherwise false.
SpecialFiberNotHyperelliptic := function(N, W)
  for p in PrimeDivisors(N) do
    M := N div p;
    if M notin GENUS0 then continue; end if;            // need X_0(M) = P^1
    Wpp := {w div p^Valuation(w,p) : w in W};            // p-coprime image of W = component group
    if GenusShimuraCurveQuotient(1, M, Wpp) ne 0 then continue; end if;
    if #Wpp gt 2 then continue; end if;                 // group quotient component: TODO
    case2 := exists{w : w in W | w mod p eq 0};
    F := GF(p^2);
    ss := ss_points(M, p, F);
    if #Wpp eq 1 then qss := ss; else qss := push_single(ss, M, Rep(Wpp diff {1}), F); end if;
    if has_compat(qss, p, F, case2) eq "no" then return true, <p, case2 select 2 else 1, M, Wpp>; end if;
  end for;
  return false, _;
end function;
