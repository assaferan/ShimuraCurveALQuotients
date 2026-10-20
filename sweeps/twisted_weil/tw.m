// Twisted Weil-polynomial test, one (D,N) level per run.  Run from the repo root:
//   MAGMA_STARTUP_FILE=/dev/null magma -b D:=6 N:=23 [IN:=sweeps/twisted_weil/all_in.txt]
//         [OUT:=sweeps/twisted_weil/out/6_23.out] sweeps/twisted_weil/tw.m < /dev/null
// (Code of the 2026-09-25 local run, scratchpad twistweil/tw.m, with only the IN/OUT defaults added.)
// Op list = twist.m's (sweeps/twisted_trace/twist.m, commit 3fda912), except that the V3 ops with 9 notin W
// (not Q-rational) are dropped entirely, and an op is kept only if M^2 = 1 on H^1(C) and it commutes with T_p
// at every tested p (otherwise it is reported NONCOMM and not tested).
// For each op h (incl. h = 1, as the reproduce-the-pipeline control) and each table prime p (g=3: p<25;
// g=4: p<=5; g=5,6: p=2), p good:  P_h(t) = P_+(t) P_-(-t), P_+- = prod (t^2 - a t + p) over T_p-eigenvalues a
// on the h = +-1 part of S_2, checked against data/hypg<g>q<p>.txt (LMFDB, hyp_count > 0).
// Sanity asserts: charpoly of T_p on K_+- is a square; P_h monic of degree 2g with p^g constant term;
// coefficient of t^(2g-1) = -Tr(T_p h | S_2) = -Tr(T M)/2 (the twisted trace of twist.m);
// P_h = P_1 mod 2 (so mod-2 / 2-rank criteria give nothing new for h != 1).
AttachSpec("ShimuraQuotients.spec");
import "ModularNonALInvolutions.m" : get_Vmu;
SetColumns(0);
if not assigned IN then IN := "sweeps/twisted_weil/all_in.txt"; end if;
if not assigned OUT then OUT := "sweeps/twisted_weil/out/" cat D cat "_" cat N cat ".out"; end if;
D := StringToInteger(D); N := StringToInteger(N); L := D*N;
tabp := func<g | g eq 3 select [2,3,5,7,11,13,17,19,23] else g eq 4 select [2,3,5] else g in {5,6} select [2] else []>;
cs := [];
for line in Split(Read(IN), "\n") do
  f := Split(line, " ");
  if #f lt 7 or line[1] eq "#" then continue; end if;
  if StringToInteger(f[3]) eq D and StringToInteger(f[4]) eq N then
    Append(~cs, <f[1], StringToInteger(f[2]), StringToInteger(f[5]), {StringToInteger(x) : x in Split(f[7], ",")}>);
  end if;
end for;
F := Open(OUT, "w");
t0 := Cputime();
// tables
TAB := AssociativeArray();
for c in cs do for p in tabp(c[3]) do
  if L mod p eq 0 or IsDefined(TAB, <c[3],p>) then continue; end if;
  TAB[<c[3],p>] := Set([l : l in Split(Read(Sprintf("data/hypg%oq%o.txt", c[3], p)), "\n") | #l gt 0]);
end for; end for;
key := func<P | "[" cat Join([IntegerToString(x) : x in Reverse(Coefficients(P))], ",") cat "]">;
MDN := ModularSymbols(L, 2, 0); SDN := CuspidalSubspace(MDN);
for p in PrimeDivisors(D) do SDN := NewSubspace(SDN, p); end for;
B := Matrix([Representation(v) : v in Basis(SDN)]); n := Nrows(B); MA := MatrixAlgebra(Rationals(), n);
chi := func<w | (-1)^#PrimeDivisors(GCD(w, D))>;
als := [Q : Q in Divisors(L) | GCD(Q, L div Q) eq 1];
ALm := AssociativeArray(); for Q in als do ALm[Q] := MA!(chi(Q)*Solution(B, B*AtkinLehnerOperator(MDN, Q))); end for;
Vfull := AssociativeArray();
if N mod 4 eq 0 then Vfull["S2"] := MA!get_Vmu(2, N, B, MDN, true); end if;
if N mod 8 eq 0 then Vfull["V2"] := MA!get_Vmu(2, N, B, MDN, false); end if;
if Valuation(N, 3) eq 2 then Vfull["V3"] := MA!get_Vmu(3, N, B, MDN, false); end if;
ps := [p : p in [2,3,5,7,11,13,17,19,23] | L mod p ne 0 and exists{c : c in cs | p in tabp(c[3])}];
Tp := AssociativeArray(); for p in ps do Tp[p] := MA!Solution(B, B*HeckeOperator(MDN, p)); end for;
fprintf F, "LEVEL %o %o dimDnew %o tall %o ncurves %o\n", D, N, n, Cputime(t0), #cs;
_<t> := PolynomialRing(Integers()); R2<x> := PolynomialRing(Rationals());
WP := function(chiT, p)   // chiT = charpoly of T on K_+- (a square c^2); returns prod (t^2 - a t + p)
  fa := Factorization(chiT); assert &and[e[2] mod 2 eq 0 : e in fa];
  c := &*[R2 | e[1]^(e[2] div 2) : e in fa];
  // P(t) = t^deg c * c(t + p/t)
  d := Degree(c); cf := Coefficients(c);
  Zt<u> := PolynomialRing(Integers());
  P := &+[Zt | Integers()!cf[i+1] * (u^2 + p)^i * u^(d-i) : i in [0..d]];
  return P;
end function;
for c in cs do
  st, id, g, W := Explode(c);
  V := VectorSpace(Rationals(), n); K := V;
  for w in W do if w ne 1 then K meet:= Kernel(ALm[w] - 1); end if; end for;
  BK := BasisMatrix(K); dK := Nrows(BK);
  if dK ne 2*g then fprintf F, "BADDIM %o %o %o\n", id, dK, g; continue; end if;
  IK := IdentityMatrix(Rationals(), dK); MK := MatrixAlgebra(Rationals(), dK);
  res := function(M) ok, S := IsConsistent(BK, BK*M); if ok then return true, MK!S; else return false, _; end if; end function;
  cand := [<"1", ALm[1], false>] cat [<"w" cat IntegerToString(Q), ALm[Q], false> : Q in als | Q notin W];
  vn := [];
  if "S2" in Keys(Vfull) and &and[IsOdd(w) : w in W] then Append(~vn, "S2"); end if;
  if "V2" in Keys(Vfull) then Append(~vn, "V2"); end if;
  if "V3" in Keys(Vfull) and 9 in W then Append(~vn, "V3"); end if;   // Q-rational V3 only
  allv := [<a, Vfull[a]> : a in vn] cat [<a cat "*" cat b, Vfull[a]*Vfull[b]> : a, b in vn | a ne b];
  for vv in allv do
    for Q in als do
      if Q in W and Q ne 1 then continue; end if;
      if "S2" in vv[1] and IsEven(Q) then continue; end if;
      Append(~cand, <vv[1] cat "*w" cat IntegerToString(Q), vv[2]*ALm[Q], false>);
    end for;
  end for;
  ops := []; nodesc := 0; notinv := [];
  for o in cand do
    ok, M := res(o[2]);
    if not ok then nodesc +:= 1; continue; end if;
    if M^2 ne IK then Append(~notinv, o[1]); continue; end if;
    if exists{y : y in ops | y[2] eq M} then continue; end if;
    Append(~ops, <o[1], M>);
  end for;
  Tk := AssociativeArray();
  for p in ps do if p in tabp(g) then ok, T := res(Tp[p]); assert ok; Tk[p] := T; end if; end for;
  tested := [p : p in Keys(Tk)]; Sort(~tested);
  // P_1 per p
  P1 := AssociativeArray();
  for p in tested do P1[p] := WP(R2!CharacteristicPolynomial(Tk[p]), p); end for;
  lines := []; fails := []; noncomm := []; minus1 := [];
  for o in ops do
    M := o[2];
    if M eq -IK then Append(~minus1, o[1]); end if;
    if exists{p : p in tested | Tk[p]*M ne M*Tk[p]} then Append(~noncomm, o[1]); continue; end if;
    Kp := Kernel(M - IK); Km := Kernel(M + IK);
    for p in tested do
      T := Tk[p];
      if Dimension(Kp) gt 0 then
        Bp := BasisMatrix(Kp); Tp_ := Solution(Bp, Bp*T); Pp := WP(R2!CharacteristicPolynomial(Tp_), p);
      else Pp := t^0; end if;
      if Dimension(Km) gt 0 then
        Bm := BasisMatrix(Km); Tm_ := Solution(Bm, Bm*T); Pm := WP(R2!CharacteristicPolynomial(Tm_), p);
      else Pm := t^0; end if;
      Ph := Pp * Evaluate(Pm, -t);
      assert Degree(Ph) eq 2*g and LeadingCoefficient(Ph) eq 1 and Coefficient(Ph, 0) eq p^g;
      assert Coefficient(Ph, 2*g-1) eq -Trace(T*M)/2;
      assert ChangeRing(Ph, GF(2)) eq ChangeRing(P1[p], GF(2));
      if o[1] eq "1" then assert Ph eq P1[p]; end if;
      inT := key(Ph) in TAB[<g,p>];
      if not inT then Append(~fails, Sprintf("%o:%o:%o", o[1], p, key(Ph))); end if;
    end for;
  end for;
  fprintf F, "RES %o %o %o %o %o %o nops %o ops %o notinv %o noncomm %o minus1 %o ps %o nfail %o fails %o\n", st, id, D, N, g,
    Join([IntegerToString(w) : w in Sort(SetToSequence(W))], ","), #ops, Join([o[1] : o in ops], ","),
    #notinv eq 0 select "-" else Join(notinv, ","), #noncomm eq 0 select "-" else Join(noncomm, ","),
    #minus1 eq 0 select "-" else Join(minus1, ","), #tested eq 0 select "-" else Join([IntegerToString(p) : p in tested], ","),
    #fails, #fails eq 0 select "-" else Join(fails, ";");
  Flush(F);
end for;
fprintf F, "DONE %o %o %o\n", D, N, Cputime(t0);
delete F;
quit;
