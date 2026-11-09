// THE SHARP DUAL TEST FOR 2 || N, D ODD, VIA THE LEVEL N/2 CURVE.
// X* = X_0^D(N)*, X'* = X_0^D(N/2)*, same scalar level M = 2DN.  With T = p_* q^* the correspondence
// X* <- Y -> X'*, the fibre count of X_0^D(N) -> X_0^D(N/2) gives, for the Heegner series Z(n) of X*
// and the plus-type series H(n) = sum_{s : -n/s^2 a discriminant} P'_{-n/s^2} of X'*,
//     T Z(n) = H(4n) + 2 H(n),
// and T is injective on J(X*)_Q when J(X*) is f-isotypic with f 2-old.  H_phi' lies in the Kohnen
// plus space S^+_{3/2}(Gamma_0(M)) (Borcherds' modularity of the vector-valued Heegner series,
// restricted to the Z/2 part of the discriminant form), and it vanishes at every index n where H(n)
// is the zero divisor.  So every Heegner component of X* lies in V = (U(4) + 2) W+, with
//     W+ = { g in S_{3/2}(Gamma_0(M)) : a_g(n) = 0 unless n = 0,3 mod 4 and some -n/s^2 carries CM on X'* },
// and a degree-zero CM divisor Delta on X* whose coefficient vector is orthogonal to V is torsion.
// Run from the repository root; the basis of S_{3/2}(Gamma_0(2436)) to precision 4002 (25 minutes) is
// cached under data/dual/ and not committed:
//   magma -b DD:=21 NN:=58 MM:=2436 BB:=4000 TGTS:="[[<-16,1>,<-120,1>,<-88,-1>,<-100,-1>]]" tools/proof-21-58/stepH-plus.m
SetColumns(0);
D := StringToInteger(DD); N := StringToInteger(NN); M := StringToInteger(MM); B := StringToInteger(BB);
error if N mod 2 ne 0 or D mod 2 eq 0, "needs D odd and N even";
Np := N div 2;
targets := eval TGTS;
function HasCM(d, D, N)
    d0 := FundamentalDiscriminant(d); _, f := IsSquare(d div d0);
    if exists{p : p in PrimeDivisors(D) | f mod p eq 0 or KroneckerSymbol(d0, p) eq 1} then return false; end if;
    if exists{p : p in PrimeDivisors(N) | f mod p ne 0 and KroneckerSymbol(d0, p) eq -1} then return false; end if;
    return true;
end function;
function Stab(d, D, N)
    e := d eq -3 select 3 else (d eq -4 select 2 else 1);
    if exists{m : m in Divisors(D*N) | m gt 1 and GCD(m, (D*N) div m) eq 1 and (d eq -4*m or (m mod 4 eq 3 and d eq -m))} then e *:= 2; end if;
    return e;
end function;
discs := func< n | [ -4*n div (r^2) : r in [1..Isqrt(4*n)] | (4*n) mod (r^2) eq 0 and ((-4*n div (r^2)) mod 4) in {0, 1} ] >;
// discriminants -n/s^2 (the plus-type index set of X'*)
pdiscs := func< n | [ -n div (s^2) : s in [1..Isqrt(n)] | n mod (s^2) eq 0 and ((-n div (s^2)) mod 4) in {0, 1} ] >;
t0 := Realtime();
cache := Sprintf("data/dual/basis_%o_%o.txt", M, B);
if System("test -f " cat cache) eq 0 then
    A0 := eval Read(cache);
else
    S := CuspidalSubspace(HalfIntegralWeightForms(Gamma0(M), 3/2));
    Bq := Basis(S, B + 2);
    A0 := Matrix(Rationals(), #Bq, B, [ Coefficient(g, n) : n in [1..B], g in Bq ]);
    Write(cache, Sprint(A0, "Magma") : Overwrite := true);
end if;
dim := Nrows(A0);
printf "M = %o: %o cusp forms, coefficients to %o from %o\n", M, dim, B, cache;
// indices where H(n) = 0 on X'*: n = 1,2 mod 4, or no -n/s^2 carries CM on X_0^D(N/2)
Hzero := [ n : n in [1..B] | (n mod 4) in {1, 2} or forall{d : d in pdiscs(n) | not HasCM(d, D, Np)} ];
printf "H vanishes at %o of %o indices\n", #Hzero, B;
AE := Matrix(Rationals(), dim, #Hzero, [ A0[i, n] : n in Hzero, i in [1..dim] ]);
K := Kernel(AE);
Wp := [ Vector(k) * A0 : k in Basis(K) ];
printf "dim W+_B = %o\n", #Wp;
// V = (U(4) + 2) W+, to precision B/4
B4 := B div 4;
V := [ [ g[4*n] + 2*g[n] : n in [1..B4] ] : g in Wp ];
// the Heegner series of X* itself vanishes at the indices with no CM point on X*: the sharp test's W;
// V meet W would be sharper still, but V alone is what the lemma gives, so use V
cmd := AssociativeArray();
for n in [1..B4] do cmd[n] := [ d : d in discs(n) | HasCM(d, D, N) ]; end for;
full := [ n : n in [1..B4] | not IsEmpty(cmd[n]) ];
cmset := {}; for n in full do for d in cmd[n] do Include(~cmset, d); end for; end for;
cmds := Sort(Setseq(cmset));
printf "%o CM discriminants of X* reachable from poles up to %o\n", #cmds, B4;
for target in targets do
    tg := AssociativeArray(); for t in target do tg[t[1]] := t[2]; end for;
    if not forall{t : t in target | t[1] in cmds} then printf "target %o: a discriminant is out of range\n", target; continue; end if;
    nrows := #V + #cmds;
    A := ZeroMatrix(Rationals(), nrows, B4); rhs := ZeroMatrix(Rationals(), nrows, 1);
    for i->v in V do for n in [1..B4] do A[i, n] := v[n]; end for; end for;
    for j->d in cmds do
        for n in full do if d in cmd[n] then A[#V + j, n] +:= 1; end if; end for;
        rhs[#V + j, 1] := IsDefined(tg, d) select Stab(d, D, N) * tg[d] else 0;
    end for;
    ok := IsConsistent(Transpose(A), Transpose(rhs));
    printf "target %o: TORSION-CERTIFIED %o (%o equations, rank %o)\n", target, ok, nrows, Rank(A);
end for;
printf "done in %o s\n", Realtime() - t0;
quit;
