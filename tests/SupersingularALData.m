// Validates the Brandt-module supersingular layer (SupersingularALData) against the
// independent D=6 geometric layer of special_fiber_modular.m: for X_0(6,1) the supersingular
// points are known explicitly (Elkies hypergeometric conic points) and the Atkin-Lehner
// involutions act as coordinate sign-flips, so the geometric fixed-point counts of w_2, w_3,
// w_6 and of Frobenius give ground truth.  The Brandt layer must reproduce all of them.
// A D=10 smoke check then exercises the intrinsic's built-in self-checks (count = genus+1,
// involution/commutation relations) on a curve with no hypergeometric model.

// --- D=6 geometric helpers (mirrors special_fiber_modular.m, which keeps them file-local) ---
d6_inert := func<p,d | not IsSplit(p, MaximalOrder(QuadraticField(d)))>;
function d6_hypg_test(p, F)
  if   p mod 24 in [19,23] then A:=19/24; B:=23/24; C:=3/2;
  elif p mod 24 in [1,5]   then A:=1/24;  B:=5/24;  C:=1/2;
  elif p mod 24 in [7,11]  then A:=7/24;  B:=11/24; C:=1/2;
  else                          A:=13/24; B:=17/24; C:=3/2; end if;
  R := PolynomialRing(F);
  return &+[ &*[F | (A+l)*(B+l)/(C+l) : l in [0..j-1]]/Factorial(j) * R.1^j : j in [0..(p div 24)]];
end function;
eqP3_test := func<a,b | a[1]*b[2]-a[2]*b[1] eq 0 and a[1]*b[3]-a[3]*b[1] eq 0 and a[2]*b[3]-a[3]*b[2] eq 0>;
function d6_ss_conic_test(p, F)
  pts := []; add := procedure(~pts,P) if not exists{q:q in pts|eqP3_test(P,q)} then Append(~pts,P); end if; end procedure;
  for r in Roots(d6_hypg_test(p,F)) do tau0:=r[1]; rz:=Sqrt(-tau0); ry:=Sqrt(-(1-tau0)/3);
    for sy in [ry,-ry] do for sz in [rz,-rz] do add(~pts,[F!1,sy,sz]); end for; end for; end for;
  if d6_inert(p,-3)  then s:=Sqrt(F!-3); add(~pts,[F!0,F!1,s]); add(~pts,[F!0,F!1,-s]); end if;
  if d6_inert(p,-4)  then s:=Sqrt(F!-1); add(~pts,[s,F!0,F!1]); add(~pts,[-s,F!0,F!1]); end if;
  if d6_inert(p,-24) then s:=Sqrt(F!-3); add(~pts,[s,F!1,F!0]); add(~pts,[-s,F!1,F!0]); end if;
  return pts;
end function;

procedure test_SupersingularALData_vs_D6_geometry()
    // geometric sign-flip Atkin-Lehner involutions on the conic x^2+3y^2+z^2:
    flip := AssociativeArray(); flip[2]:=[1,-1,1]; flip[3]:=[-1,1,1]; flip[6]:=[-1,-1,1];
    for p in [q : q in PrimesInInterval(7, 59) | 6 mod q ne 0] do
        F := GF(p^2);
        ss := d6_ss_conic_test(p, F);
        geomfix := func<e | #[Q : Q in ss | eqP3_test([e[1]*Q[1],e[2]*Q[2],e[3]*Q[3]], Q)]>;
        geom_frob := #[Q : Q in ss | eqP3_test([Frobenius(Q[1]),Frobenius(Q[2]),Frobenius(Q[3])], Q)];

        n, frob, al := SupersingularALData(6, p);
        assert n eq #ss;                                  // same supersingular count
        assert #Fix(frob) eq geom_frob;                   // Frobenius = w_p (F_p-rational count)
        for q in [2,3,6] do
            assert #Fix(al[q]) eq geomfix(flip[q]);        // Atkin-Lehner fixed-point counts
        end for;
    end for;
end procedure;

procedure test_SupersingularALData_D10_smoke()
    // No geometric oracle here; rely on the intrinsic's internal self-checks (count = genus+1,
    // involution/commutation asserts) plus a couple of externally-checkable facts.
    for p in [q : q in PrimesInInterval(7, 47) | 10 mod q ne 0] do
        n, frob, al := SupersingularALData(10, p);
        assert n eq GenusShimuraCurveQuotient(10, p, {Integers()|1}) + 1;
        assert Keys(al) eq {2, 5, 10};
        assert al[10] eq al[2]*al[5];                     // w_10 = w_2 w_5
        assert IsEven(n - #Fix(frob));                    // non-rational ss points pair up under Frobenius
    end for;
end procedure;

test_SupersingularALData_vs_D6_geometry();
test_SupersingularALData_D10_smoke();
