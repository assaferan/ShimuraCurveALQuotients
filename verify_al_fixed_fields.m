// verify_al_fixed_fields.m
//
// Decide between candidate models of an Atkin-Lehner quotient C = X_0(D,N)/W when point counts
// cannot: two curves with isogenous Jacobians have the SAME a_p at every good prime, so the
// trace formula pins only the isogeny class.  The fixed points of the AL involutions, and their
// fields of definition, are invariants of the CURVE (not of its Jacobian), and do decide.
//
// For an even sextic model y^2 = A x^6 + B x^4 + C x^2 + D the residual group W_full/W acts by
//     sigma      : (x,y) -> (-x,  y)   Fix = (0, +-sqrt(D))            field Q(sqrt(D))
//     sigma*iota : (x,y) -> (-x, -y)   Fix = 2 pts at infty, y/x^3 = +-sqrt(A)   field Q(sqrt(A))
//     iota       : (x,y) -> ( x, -y)   Fix = the 6 Weierstrass points (roots of the sextic)
// The square classes of A and D are invariants of the model up to x -> lambda x, y -> nu y.
//
// Each model involution is matched to its AL class by comparing the a_p of its quotient with the
// trace formula for X/(W join class) -- no hand assignment.
//
// Usage:   magma -b verify_al_fixed_fields.m
// Edit the PARAMETERS block to point at a different quotient or different candidates.

SetColumns(0);
AttachSpec("ShimuraQuotients.spec");
P<x> := PolynomialRing(Rationals());

// ---------------------------------------------------------------- PARAMETERS
D := 34;
N := 3;
Wtarget := {1, 102};                       // the quotient whose model we are testing
primes := [5,7,11,13,19,23,29,31,41,43,47];   // good primes used for the a_p matching
cands := [
  <"cand1", -1152/37*x^6 + 2709/37*x^4 + 2529/37*x^2 - 864/37>,
  <"cand2",  36*x^6 + 117*x^4 + 18*x^2 - 27>
];
// ---------------------------------------------------------------------------

function sqfree(q)
  n := Numerator(q)*Denominator(q); s := Sign(n); n := AbsoluteValue(n); r := 1;
  for t in Factorization(n) do if IsOdd(t[2]) then r *:= t[1]; end if; end for;
  return s*r;
end function;

function fldname(F)
  if Degree(F) eq 1 then return "Q"; end if;
  if Degree(F) eq 2 then return Sprintf("Q(sqrt(%o))", sqfree(Discriminant(MaximalOrder(F)))); end if;
  return Sprintf("deg %o, disc %o", Degree(F), Factorization(Discriminant(MaximalOrder(F))));
end function;

// y^2 = A u^3 + B u^2 + C u + D  ==>  Y^2 = U^3 + B U^2 + (A C) U + A^2 D,  U = A u, Y = A y
function ecfromcubic(A,B,C,D)
  return MinimalModel(EllipticCurve([0, B, 0, A*C, A^2*D]));
end function;

curves := GetHyperellipticCandidates();
C := rep{Y : Y in curves | Y`D eq D and Y`N eq N and Y`W eq Wtarget};
DN := D*N;
Wfull := Set(Divisors(DN));

printf "C = X_0(%o,%o)/%o,  genus %o\n\n", D, N, Sort([w : w in Wtarget]), C`g;

// ---- the AL classes acting on C: cosets of Wtarget in Wfull -----------------
classes := [];
seen := {};
for w in Sort([t : t in Wfull]) do
  if w in seen then continue; end if;
  cl := {AtkinLehnerMul(w, m, DN) : m in Wtarget};
  seen join:= cl;
  if cl ne Wtarget then Append(~classes, Sort([t : t in cl])); end if;
end for;

// ---- EXPECTED table --------------------------------------------------------
printf "================ EXPECTED (from the Shimura curve) ================\n";
printf "%-16o %-6o %-26o %o\n", "AL class", "#Fix", "CM disc (pts upstairs)", "field(s) of definition on C";
expected := AssociativeArray();      // class -> list of <disc, [fields]>
quotgenus := AssociativeArray();
for cl in classes do
  nfix := CountFixedPointsOnQuotient(cl[1], C);
  Wq := Wtarget join Set(cl);
  Xq := rep{Y : Y in curves | Y`D eq D and Y`N eq N and Y`W eq Wq};
  quotgenus[cl] := Xq`g;
  rows := [* *];
  first := true;
  for m in cl do
    e := NumFixedPointsByCMOrder(D, N, m);
    for d in Sort([k : k in Keys(e)]) do
      if e[d] eq 0 then continue; end if;
      fs := FieldsOfDefinitionOfCMPoint(C, d);
      Append(~rows, <d, fs>);
      printf "%-16o %-6o %-26o %o\n",
             first select Sprintf("{w_%o,w_%o}", cl[1], cl[2]) else "",
             first select nfix else "",
             Sprintf("%o  (%o from w_%o)", d, e[d], m),
             [fldname(F) : F in fs];
      first := false;
    end for;
  end for;
  expected[cl] := rows;
  printf "%-16o %-6o quotient X/%o has genus %o\n\n", "", "", Sort([t : t in Wq]), Xq`g;
end for;

// ---- a_p of the target and of each genus-1 quotient, for matching ----------
tgt_ap := [p+1-ComputePointsViaTrace(C,p,1) : p in primes | DN mod p ne 0];
qt_ap := AssociativeArray();
for cl in classes do
  if quotgenus[cl] ne 1 then continue; end if;
  Wq := Wtarget join Set(cl);
  Xq := rep{Y : Y in curves | Y`D eq D and Y`N eq N and Y`W eq Wq};
  qt_ap[cl] := [p+1-ComputePointsViaTrace(Xq,p,1) : p in primes | DN mod p ne 0];
end for;

// ---- OBSERVED table, per candidate ----------------------------------------
for cd in cands do
  name := cd[1]; f := cd[2];
  A := Coefficient(f,6); B := Coefficient(f,4); Cc := Coefficient(f,2); Dd := Coefficient(f,0);
  Cv := HyperellipticCurve(f);
  printf "================ OBSERVED: %o ================\n", name;
  printf "y^2 = %o\n", f;
  require Genus(Cv) eq C`g : "candidate has the wrong genus";

  ap := [p+1-#Points(ChangeRing(Cv,GF(p))) : p in primes | DN mod p ne 0];
  printf "a_p matches the trace formula: %o   (so point counts alone cannot decide)\n\n", ap eq tgt_ap;

  // match each model involution to an AL class via the a_p of its quotient
  Ea := ecfromcubic(A,B,Cc,Dd);      // sigma      : (x,y) -> (-x, y)
  Eb := ecfromcubic(Dd,Cc,B,A);      // sigma*iota : (x,y) -> (-x,-y)
  invs := [* <"sigma      (x,y)->(-x,y) ", Ea, Dd>,
             <"sigma*iota (x,y)->(-x,-y)", Eb, A> *];

  npass := 0; nfail := 0;
  for iv in invs do
    lbl := iv[1]; E := iv[2]; coef := iv[3];
    aE := [TraceOfFrobenius(E,p) : p in primes | DN mod p ne 0];
    cl := 0;
    for c2 in classes do
      if quotgenus[c2] eq 1 and qt_ap[c2] eq aE then cl := c2; end if;
    end for;
    if cl cmpeq 0 then
      printf "  %o -> quotient (cond %o) matches NO AL class\n", lbl, Conductor(E); nfail +:= 1; continue;
    end if;
    obs := QuadraticField(sqfree(coef));
    exp := [r : r in expected[cl]];
    okrow := &or[ &or[ IsIsomorphic(obs, F) : F in r[2] ] : r in exp ];
    printf "  %o -> AL class {w_%o,w_%o} (quotient cond %o)\n", lbl, cl[1], cl[2], Conductor(E);
    printf "        Fix field observed  : %o\n",
           (sqfree(coef) eq 1) select "Q" else Sprintf("Q(sqrt(%o))", sqfree(coef));
    printf "        Fix field expected  : %o\n", &cat[Sprintf("%o ", [fldname(F) : F in r[2]]) : r in exp];
    printf "        %o\n", okrow select "MATCH" else "*** CONTRADICTION ***";
    if okrow then npass +:= 1; else nfail +:= 1; end if;
  end for;

  // the hyperelliptic class: Weierstrass points = roots of the sextic
  hcl := 0;
  for c2 in classes do if quotgenus[c2] eq 0 then hcl := c2; end if; end for;
  printf "  iota       (x,y)->(x,-y) -> AL class {w_%o,w_%o} (hyperelliptic, quotient genus 0)\n", hcl[1], hcl[2];
  expfields := &cat[[F : F in r[2]] : r in expected[hcl]];
  for t in Factorization(f) do
    if Degree(t[1]) eq 0 then continue; end if;
    F := NumberField(t[1] : DoLinearExtension);
    ok := &or[ IsIsomorphic(F, G) : G in expfields ];
    printf "        factor %-26o -> %-28o %o\n", t[1], fldname(F),
           ok select "MATCH" else "*** CONTRADICTION ***";
    if ok then npass +:= 1; else nfail +:= 1; end if;
  end for;
  printf "\n  VERDICT for %o:  %o rows match, %o contradict\n\n", name, npass, nfail;
end for;
quit;
