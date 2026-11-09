// THE 2 || N TEST FOR DIVISORS IN ODD DISCRIMINANTS, THROUGH THE LEVEL N/2 CURVE.
// A degree-zero CM divisor Delta on X* = X_0^D(N)* is pushed to X'* = X_0^D(N/2)* by the correspondence
// T of the lemma (T C(d) = (1 + chi_d(2)) C'(d) + (2 - chi_{d/4}(2)) C'(d/4), orbifold divisors C = P/e),
// written as a combination of the plus-type Heegner divisors H(n) = sum_s C'(-n/s^2) of X'* (in this
// dictionary every discriminant has its own index n = |d|), and tested against W+ (the plus-space cusp
// forms vanishing where H vanishes), which contains every component H_phi'.  T is injective on
// Jac(X*)_Q when the Jacobian is 2-old, so Delta is torsion iff T Delta is.
// Run from the repository root, with the cache written by stepH-plus.m:
//   magma -b DD:=39 NN:=34 MM:=2652 BB:=4000 CACHE:=data/dual/basis_2652_4000.txt TGTS:="[[<-15,1>,<-52,1>,<-84,-1>,<-132,-1>]]" tools/proof-39-34/stepH-odd.m
// (targets list the multiplicities of the POINTS P_d of X*, as everywhere else)
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
// chi_d(2): the Kronecker symbol (d/2) if 2 does not divide the conductor, 0 if it does
function Chi2(d)
    d0 := FundamentalDiscriminant(d); _, f := IsSquare(d div d0);
    return f mod 2 eq 0 select 0 else KroneckerSymbol(d0, 2);
end function;
pdiscs := func< n | [ -n div (s^2) : s in [1..Isqrt(n)] | n mod (s^2) eq 0 and ((-n div (s^2)) mod 4) in {0, 1} ] >;
t0 := Realtime();
A0 := eval Read(CACHE);
dim := Nrows(A0);
printf "M = %o: %o cusp forms, coefficients to %o from %o\n", M, dim, B, CACHE;
Hzero := [ n : n in [1..B] | (n mod 4) in {1, 2} or forall{d : d in pdiscs(n) | not HasCM(d, D, Np)} ];
K := Kernel(Matrix(Rationals(), dim, #Hzero, [ A0[i, n] : n in Hzero, i in [1..dim] ]));
Wp := [ Vector(k) * A0 : k in Basis(K) ];
printf "dim W+_B = %o\n", #Wp;
plusidx := [ n : n in [1..B] | n notin Hzero ];                 // indices where H(n) is a nonzero divisor
hd := AssociativeArray();
for n in plusidx do hd[n] := [ d : d in pdiscs(n) | HasCM(d, D, Np) ]; end for;
dprime := Sort(Setseq(&join{ Set(hd[n]) : n in plusidx }));     // discriminants of X'* reachable
for target in targets do
    // C-units on X*: alpha_d = e_d * m_d; then T
    mu := AssociativeArray();
    ok := true;
    for t in target do
        d := t[1]; alpha := Stab(d, D, N) * t[2];
        if not HasCM(d, D, N) then printf "target %o: %o carries no CM point on X*\n", target, d; ok := false; break; end if;
        if HasCM(d, D, Np) then mu[d] := (IsDefined(mu, d) select mu[d] else 0) + alpha * (1 + Chi2(d)); end if;
        if d mod 4 eq 0 and ((d div 4) mod 4) in {0, 1} and HasCM(d div 4, D, Np) then
            mu[d div 4] := (IsDefined(mu, d div 4) select mu[d div 4] else 0) + alpha * (2 - Chi2(d div 4));
        end if;
    end for;
    if not ok then continue; end if;
    if not forall{d : d in Keys(mu) | d in dprime} then printf "target %o: pushed discriminant out of range\n", target; continue; end if;
    // realise T Delta = sum_n c'(n) H(n): one equation per reachable discriminant of X'*
    nrows := #Wp + #dprime;
    A := ZeroMatrix(Rationals(), nrows, #plusidx); rhs := ZeroMatrix(Rationals(), nrows, 1);
    for i->w in Wp do for j->n in plusidx do A[i, j] := w[n]; end for; end for;
    for k->d in dprime do
        for j->n in plusidx do if d in hd[n] then A[#Wp + k, j] +:= 1; end if; end for;
        rhs[#Wp + k, 1] := IsDefined(mu, d) select mu[d] else 0;
    end for;
    cons := IsConsistent(Transpose(A), Transpose(rhs));
    printf "target %o -> T Delta = %o: TORSION-CERTIFIED %o\n", target, [ <d, mu[d]> : d in Sort(Setseq(Keys(mu))) ], cons;
end for;
printf "done in %o s\n", Realtime() - t0;
quit;
