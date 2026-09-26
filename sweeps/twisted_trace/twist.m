// Twisted trace (twisted Lefschetz) non-hyperellipticity test, one (D,N) level per run.
// Adapted from v3cert/probe6.m; this is the code of the 2026-09-25 local counting run
// (scratchpad twistcount/twist.m) with I/O changed (paths are args, output is per-level) and one rule
// changed: V3 ops are used at every good p when 9 in W ("v3all", below).
//
// Run from the repo root:
//   magma -b D:=15 N:=146 [BOUND:=weil|pb] [PB:=59] [IN:=sweeps/twisted_trace/curves_pending.txt]
//         [OUT:=sweeps/twisted_trace/out/15_146.out] sweeps/twisted_trace/twist.m < /dev/null
//
// For each curve C = X_0^D(N)/W listed in IN at this (D,N) and each h in the known group:
//   residual ALs, and V*W_Q (V in S2,V2,V3 and pairwise products, same descent conditions as
//   CheckModularNonALInvolutionModSym), plus h = 1,
// and each q = p^v with p good (range set by BOUND below), record every |Tr(s_v(T_p) h)/2| > q+1 on H^1(C).
// (report.py counts only tr < -(q+1); tr > q+1 is impossible for any curve and is flagged as a BUG.)
// Ops involving V3: if 9 in W they are used at EVERY good p (rule "v3all"): Galois acts by
// sigma(V3) = V3*W9, so V3 (and V3*w) induces an involution of C defined over Q when 9 in W.  If 9 notin W,
// V3 need not be Q-rational and those ops are only used at p = 1 mod 3 (report.py never counts them).
// (Before v3all, V3 ops were used at p = 1 mod 3 only, for every W; RES lines of that rule carry no v3all.)
// Any op not commuting with T_p at p is skipped (counted in noncomm): the safety net for every op.
// V3ONLY:=1 (supplement mode): only curves with 9 in W, only the V3-containing ops, only p = 2 mod 3 --
// exactly what a pre-v3all run left untested.  Default OUT suffix .v3.  OUT is truncated at start; a level
// is complete iff OUT ends with a DONE line.
// Output: LEVEL / RES (one per curve) / BADDIM / DONE lines.
AttachSpec("ShimuraQuotients.spec");
import "ModularNonALInvolutions.m" : get_Vmu;
SetColumns(0);
// Prime-power range.  BOUND:=weil (default): per curve of genus g, every q = p^v < 4g^2 (optionally
// also p <= PB if PB is given).  This is exhaustive: h is a Q-rational involution, so it commutes with
// Frob and every eigenvalue of Frob_q h on H^1 has |.| = sqrt(q); hence |tr| <= 2g sqrt(q), and
// tr < -(q+1) forces q+1 < 2g sqrt(q), i.e. q < (g + sqrt(g^2-1))^2 < 4g^2.  (Same bound as FilterByTrace.)
// BOUND:=pb: the 2026-09-25 local run exactly: p <= PB (default 59), q <= PB^2, for every curve.
if not assigned BOUND then BOUND := "weil"; end if;
error if BOUND notin {"weil", "pb"}, "BOUND must be weil or pb";
if not assigned PB then PB := BOUND eq "pb" select "59" else "0"; end if;   // 0 = no prime cap
// PMIN (supplement mode): only primes p >= PMIN.  Used to fill the q < 4g^2 gaps left by a PB=59 run
// (report.py writes levels_supplement.txt with the PMIN each level needs).  Default 0 = no lower cap.
if not assigned PMIN then PMIN := "0"; end if;
if not assigned V3ONLY then V3ONLY := "0"; end if;
error if V3ONLY notin {"0", "1"}, "V3ONLY must be 0 or 1";
V3ONLY := V3ONLY eq "1";
if not assigned IN then IN := "sweeps/twisted_trace/curves_pending.txt"; end if;
if not assigned OUT then OUT := "sweeps/twisted_trace/out/" cat D cat "_" cat N cat (V3ONLY select ".v3" else (PMIN eq "0" select "" else ".supp")) cat ".out"; end if;
PMIN := StringToInteger(PMIN);
D := StringToInteger(D); N := StringToInteger(N); PB := StringToInteger(PB);
L := D*N;
cs := [];
for line in Split(Read(IN), "\n") do
  f := Split(line, " ");
  if #f lt 7 or line[1] eq "#" then continue; end if;
  if StringToInteger(f[3]) eq D and StringToInteger(f[4]) eq N and (not V3ONLY or 9 in {StringToInteger(x) : x in Split(f[7], ",")}) then
    Append(~cs, <f[1], StringToInteger(f[2]), StringToInteger(f[5]), {StringToInteger(x) : x in Split(f[7], ",")}>);
  end if;
end for;
// qmax(g): largest q tested on a curve of genus g.
qmax := func<g | BOUND eq "pb" select PB^2 else 4*g^2 - 1>;
pcap := #cs eq 0 select 1 else Max([qmax(c[3]) : c in cs]);
if PB gt 0 then pcap := Min(pcap, PB); end if;
t0 := Cputime();
MDN := ModularSymbols(L, 2, 0); SDN := CuspidalSubspace(MDN);
for p in PrimeDivisors(D) do SDN := NewSubspace(SDN, p); end for;
B := Matrix([Representation(v) : v in Basis(SDN)]); n := Nrows(B); MA := MatrixAlgebra(Rationals(), n);
tsp := Cputime(t0);
chi := func<w | (-1)^#PrimeDivisors(GCD(w, D))>;
als := [Q : Q in Divisors(L) | GCD(Q, L div Q) eq 1];
ALm := AssociativeArray(); for Q in als do ALm[Q] := MA!(chi(Q)*Solution(B, B*AtkinLehnerOperator(MDN, Q))); end for;
Vfull := AssociativeArray();
if N mod 4 eq 0 then Vfull["S2"] := MA!get_Vmu(2, N, B, MDN, true); end if;
if N mod 8 eq 0 then Vfull["V2"] := MA!get_Vmu(2, N, B, MDN, false); end if;
if Valuation(N, 3) eq 2 then Vfull["V3"] := MA!get_Vmu(3, N, B, MDN, false); end if;
ps := [p : p in PrimesUpTo(pcap) | L mod p ne 0 and p ge PMIN and (not V3ONLY or p mod 3 eq 2)];
Tp := AssociativeArray(); for p in ps do Tp[p] := MA!Solution(B, B*HeckeOperator(MDN, p)); end for;
tall := Cputime(t0);
F := Open(OUT, "w");
fprintf F, "LEVEL %o %o dimDnew %o tspace %o tall %o ncurves %o bound %o PB %o pmax %o pmin %o v3all 1 v3only %o\n", D, N, n, tsp, tall, #cs, BOUND, PB, pcap, PMIN, V3ONLY select 1 else 0;
for c in cs do
  st, id, g, W := Explode(c);
  V := VectorSpace(Rationals(), n); K := V;
  for w in W do if w ne 1 then K meet:= Kernel(ALm[w] - 1); end if; end for;
  BK := BasisMatrix(K); dK := Nrows(BK);
  if dK ne 2*g then fprintf F, "BADDIM %o %o %o\n", id, dK, g; continue; end if;
  IK := IdentityMatrix(Rationals(), dK);
  res := function(M) ok, S := IsConsistent(BK, BK*M); if ok then return true, MatrixAlgebra(Rationals(), dK)!S; else return false, _; end if; end function;
  // candidate ops: <name, matrix, restricted to p = 1 mod 3 (V3 op with 9 notin W)>
  cand := [<"1", ALm[1], false>] cat [<"w" cat IntegerToString(Q), ALm[Q], false> : Q in als | Q notin W];
  vn := [];
  if "S2" in Keys(Vfull) and &and[IsOdd(w) : w in W] then Append(~vn, "S2"); end if;
  if "V2" in Keys(Vfull) then Append(~vn, "V2"); end if;
  if "V3" in Keys(Vfull) then
    nc := false;
    if 9 notin W then nc := exists{w : w in W | (w div 3^Valuation(w, 3)) mod 3 eq 2}; end if;
    if not nc then Append(~vn, "V3"); end if;
  end if;
  allv := [<a, Vfull[a]> : a in vn] cat [<a cat "*" cat b, Vfull[a]*Vfull[b]> : a, b in vn | a ne b];
  for vv in allv do
    for Q in als do
      if Q in W and Q ne 1 then continue; end if;
      if "S2" in vv[1] and IsEven(Q) then continue; end if;
      if "V3" in vv[1] and 9 notin W and exists{p : p in PrimeDivisors(Q) | (p^Valuation(Q,p) mod 3) eq 2} then continue; end if;
      Append(~cand, <vv[1] cat "*w" cat IntegerToString(Q), vv[2]*ALm[Q], "V3" in vv[1] and 9 notin W>);
    end for;
  end for;
  ops := []; nodesc := 0;
  for o in cand do
    ok, M := res(o[2]);
    if not ok then nodesc +:= 1; continue; end if;
    if exists{x : x in ops | x[2] eq M} then continue; end if;
    Append(~ops, <o[1], M, o[3]>);
  end for;
  if V3ONLY then ops := [o : o in ops | "V3" in o[1]]; end if;
  viol := []; noncomm := 0; h1viol := false;
  QM := qmax(g);
  for p in ps do
    if p gt QM then continue; end if;
    ok, T := res(Tp[p]); assert ok;
    s0 := 2*IK; s1 := T; v := 1;
    comm := [T*o[2] eq o[2]*T : o in ops];
    noncomm +:= #[x : x in comm | not x];
    while p^v le QM do
      q := p^v;
      for i->o in ops do
        if not comm[i] then continue; end if;
        if o[3] and p mod 3 ne 1 then continue; end if;
        tr := Trace(s1*o[2]) / 2;
        if Abs(tr) gt q + 1 then Append(~viol, <q, o[1], tr>); if o[1] eq "1" then h1viol := true; end if; end if;
      end for;
      s2 := T*s1 - p*s0; s0 := s1; s1 := s2; v +:= 1;
    end while;
  end for;
  fprintf F, "RES %o %o %o %o %o %o nops %o nodesc %o noncomm %o h1 %o nviol %o viol %o ops %o\n", st, id, D, N, g,
     Join([IntegerToString(w) : w in Sort(SetToSequence(W))], ","), #ops, nodesc, noncomm, h1viol, #viol,
     (#viol eq 0 select "-" else &cat[Sprintf("%o:%o:%o;", x[1], x[2], x[3]) : x in viol]), (#ops eq 0 select "-" else Join([o[1] : o in ops], ",")) cat " qmax " cat IntegerToString(QM) cat " pmax " cat IntegerToString(Max([0] cat [p : p in ps | p le QM])) cat " pmin " cat IntegerToString(PMIN) cat " v3all 1 v3only " cat (V3ONLY select "1" else "0");
  Flush(F);
end for;
fprintf F, "DONE %o %o %o\n", D, N, Cputime(t0);
delete F;
quit;
