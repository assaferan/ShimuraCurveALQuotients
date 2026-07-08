// Independent hyperellipticity oracle for AL quotients of classical modular curves (D = 1).
//
// For X_0(N)/W of genus g >= 3, the canonical model lives in P^{g-1} via the W-invariant
// weight-2 cusp forms.  By Petri's theorem the number of independent quadrics through the
// canonical image is
//        (g-2)(g-3)/2   iff the curve is NON-hyperelliptic,
//        (g-1)(g-2)/2   iff the curve is hyperelliptic
// (the latter image being a rational normal curve).  Counting quadric relations among the
// q-expansions therefore decides hyperellipticity outright, independently of any trace /
// Weil / special-fiber argument -- so it serves as a ground-truth cross-check for the D = 1
// part of the pipeline.  The check is rigorous: relations among the (weight-4) products are
// tested up to the weight-4 Sturm bound for Gamma_0(N).
//
// Calibrated against known answers (X_0(42)/<w_2>, X_0(42)/<w_7> non-hyp; X_0(30), X_0(40),
// X_0(47) hyp) and run across all D = 1, g >= 3 curves: it agrees with all 59 pipeline
// hyperelliptic determinations and all 14 special-fiber non-hyperelliptic determinations,
// and resolved the 43 then-undetermined curves (40 non-hyp, 3 hyp).
//
// Only valid for D = 1 (needs q-expansions); Shimura curves D > 1 require the special-fiber
// method (special_fiber_modular.m).

// Returns "HYP" / "NON-HYP" / "low_genus" / "??", the quotient genus, and the quadric count.
D1QuotientHyperellipticType := function(N, W)
  S := CuspForms(N);
  g := Dimension(S);
  prec := (4*Index(Gamma0(N)) div 12) + 12;       // weight-4 Sturm bound + buffer (rigorous)
  bas := [qExpansion(f, prec) : f in Basis(S)];
  inv := VectorSpace(Rationals(), g);
  for q in W do
    if q eq 1 then continue; end if;
    inv := inv meet Eigenspace(Matrix(Rationals(), AtkinLehnerOperator(S, q)), 1);
  end for;
  gw := Dimension(inv);
  if gw lt 3 then return "low_genus", gw, 0; end if;
  forms := [ &+[B[i]*bas[i] : i in [1..g]] : B in Basis(inv) ];
  prods := [ forms[i]*forms[j] : i in [1..gw], j in [1..gw] | i le j ];
  cols := prec - 2;
  Mt := Matrix(Rationals(), #prods, cols, [[Coefficient(p,k) : k in [1..cols]] : p in prods]);
  nq := #prods - Rank(Mt);
  nonhyp := (gw-2)*(gw-3) div 2; hyp := (gw-1)*(gw-2) div 2;
  return (nq eq nonhyp select "NON-HYP" else (nq eq hyp select "HYP" else "??")), gw, nq;
end function;
