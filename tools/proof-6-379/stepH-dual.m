// THE DUAL CERTIFICATE OF X_0^6(379)*.  Does a weakly holomorphic form f of weight 1/2 on Gamma_0(M)
// with the theta multiplier, holomorphic away from the cusp oo, exist whose Borcherds form on
// X_0^D(N)* has a prescribed divisor?  Unknowns: the principal part c(-n), 1 <= n <= B.  Constraints:
//   (i)  f exists iff its principal part pairs to zero with every weight-3/2 cusp form g on
//        Gamma_0(M) (theta^3 multiplier): the residue theorem for f g dtau on X_0(M);
//   (ii) the divisor is sum_n c(-n) sum_r (1/e_d) P_d over d = -4n/r^2 with CM(d) nonempty
//        (Guo-Yang, Lemma 25), so for every such d up to 4B the inner sum must equal e_d times the
//        prescribed multiplicity (0 for d not in the target).
// A rational solution, scaled to integers, is a certificate that a Borcherds form with that
// divisor exists; the form itself is never constructed.
// Run from the repository root (the basis of S_{3/2}(Gamma_0(4548)) to precision 2402 takes about four
// hours and is cached under data/dual/, not committed):
//   magma -b DD:=6 NN:=379 MM:=4548 BB:=2400 CAP:=700 TGT:="[<-123,1>,<-163,1>,<-52,-1>,<-168,-1>]" tools/proof-6-379/stepH-dual.m
// CAP restricts the poles to order <= CAP: at 700 the system is consistent, at 500 it is not; the
// other fibre pairs ({-75,-228}, {-43,-312}, {-147,-408} against {-52,-168}) are consistent too, and
// P_-123 - P_-163 (two points of one fibre) is refused.  The output line CERT := [...] is the
// certificate, banked primitive in data/dual/cert_6_379.m.
SetColumns(0);
D := StringToInteger(DD); N := StringToInteger(NN); M := StringToInteger(MM); B := StringToInteger(BB);
target := eval TGT;
// CAP: restrict the unknowns to poles of order <= CAP (<= B) to look for a certificate of small pole order
cap := assigned CAP select StringToInteger(CAP) else B;
// CM(d) nonempty on X_0^D(N), N squarefree: at p | D, p does not divide the conductor and p does not
// split in Q(sqrt d); at p | N, p divides the conductor or p is not inert.
function HasCM(d, D, N)
    d0 := FundamentalDiscriminant(d); _, f := IsSquare(d div d0);
    if exists{p : p in PrimeDivisors(D) | f mod p eq 0 or KroneckerSymbol(d0, p) eq 1} then return false; end if;
    if exists{p : p in PrimeDivisors(N) | f mod p ne 0 and KroneckerSymbol(d0, p) eq -1} then return false; end if;
    return true;
end function;
// e_d: the stabiliser order of a point of CM(d) in N_B^+(O)/Q^*: the elliptic order times 2 if an
// Atkin-Lehner involution fixes it (Ogg: w_m fixes discriminants -4m, and -m when m = 3 mod 4).
function Stab(d, D, N)
    e := d eq -3 select 3 else (d eq -4 select 2 else 1);
    if exists{m : m in Divisors(D*N) | m gt 1 and GCD(m, (D*N) div m) eq 1 and (d eq -4*m or (m mod 4 eq 3 and d eq -m))} then e *:= 2; end if;
    return e;
end function;
discs := func< n | [ -4*n div (r^2) : r in [1..Isqrt(4*n)] | (4*n) mod (r^2) eq 0 and ((-4*n div (r^2)) mod 4) in {0, 1} ] >;
t0 := Realtime();
// the basis is the expensive part (an hour or more); cache its coefficient matrix per (M, precision)
cache := Sprintf("data/dual/basis_%o_%o.txt", M, B);   // not committed: 7 MB, hours to recompute
if System("test -f " cat cache) eq 0 then
    A0 := eval Read(cache);
    Bas := [ [ A0[i, n] : n in [1..cap] ] : i in [1..Nrows(A0)] ];
    printf "M = %o: dim S_{3/2} = %o, coefficients to %o read from %o\n", M, #Bas, B, cache;
else
    S := CuspidalSubspace(HalfIntegralWeightForms(Gamma0(M), 3/2));
    Bq := Basis(S, B + 2);
    Bas := [ [ Coefficient(g, n) : n in [1..B] ] : g in Bq ];
    A0 := Matrix(Rationals(), #Bas, B, &cat Bas);
    SetColumns(0); Write(cache, Sprint(A0, "Magma") : Overwrite := true);
    printf "M = %o: dim S_{3/2} = %o, basis to precision %o in %o s (cached in %o)\n", M, #Bas, B + 2, Realtime() - t0, cache;
    Bas := [ b[1..cap] : b in Bas ];
end if;
B := cap;
// rows: one per cusp form (obstruction), one per CM discriminant d with |d| <= 4B (divisor)
cmset := {}; for n in [1..B] do for d in discs(n) do if HasCM(d, D, N) then Include(~cmset, d); end if; end for; end for;
cmds := Sort(Setseq(cmset));
printf "%o CM discriminants reachable from poles up to %o\n", #cmds, B;
tg := AssociativeArray(); for t in target do tg[t[1]] := t[2]; end for;
error if not forall{t : t in target | t[1] in cmds}, "a target discriminant carries no CM point or is out of range";
nrows := #Bas + #cmds;
A := ZeroMatrix(Rationals(), nrows, B); rhs := ZeroMatrix(Rationals(), nrows, 1);
for i->g in Bas do for n in [1..B] do A[i, n] := g[n]; end for; end for;
for j->d in cmds do
    for n in [1..B] do if d in discs(n) then A[#Bas + j, n] +:= 1; end if; end for;
    rhs[#Bas + j, 1] := IsDefined(tg, d) select Stab(d, D, N) * tg[d] else 0;
end for;
printf "system: %o equations, %o unknowns, rank %o (%o s)\n", nrows, B, Rank(A), Realtime() - t0;
ok, sol := IsConsistent(Transpose(A), Transpose(rhs));
printf "CONSISTENT: %o\n", ok;
if ok then
    c := Eltseq(sol);
    den := LCM([Denominator(x) : x in c]); c := [Integers() | den * x : x in c];
    supp := [n : n in [1..B] | c[n] ne 0];
    printf "a certificate with %o nonzero coefficients, pole order %o, scaled by %o; coefficients at the target indices: %o\n",
           #supp, Maximum(supp), den, [ <n, c[n]> : n in supp | not forall{d : d in discs(n) | not HasCM(d, D, N)} ];
    printf "CERT := %o;\n", c;
end if;
printf "RESULT %o\n", ok;
quit;
