// Which CM divisors on X_0^15(38)* supported on a given set of discriminants are divisors of Borcherds
// products (c orthogonal to ALL weight-3/2 cusp forms at level 1140, poles <= 4000)?  And which are torsion
// by the sharp test (c orthogonal to W_f, the cusp forms vanishing at the empty indices in the a_7 = -4
// eigenspace)?  The gap between the two is what the certificate solve for Delta ran into.
SetColumns(0);
D := 15; N := 38; M := 1140; B := 4000;
S5 := [-228, -760, -12, -40, -532];
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
A0 := eval Read("basis_1140_4000.txt");
dim := Nrows(A0);
cmd := AssociativeArray();
for n in [1..B] do cmd[n] := [ d : d in discs(n) | HasCM(d, D, N) ]; end for;
cmds := Sort(Setseq(&join{ Set(cmd[n]) : n in [1..B] }));
others := [ d : d in cmds | d notin S5 ];
printf "%o CM discriminants; stabilisers on S5: %o\n", #cmds, [ <d, Stab(d, D, N)> : d in S5 ];
// c with c ⊥ S and zero divisor outside S5: kernel of [A0 ; rows of the other discriminants]
nrows := dim + #others;
A := ZeroMatrix(Rationals(), nrows, B);
for i in [1..dim] do for n in [1..B] do A[i, n] := A0[i, n]; end for; end for;
for j->d in others do for n in [1..B] do if d in cmd[n] then A[dim + j, n] +:= 1; end if; end for; end for;
K := Kernel(Transpose(A));   // vectors c (as rows) with A c^T = 0
printf "kernel dimension %o\n", Dimension(K);
Dv := ZeroMatrix(Rationals(), B, #S5);
for j->d in S5 do for n in [1..B] do if d in cmd[n] then Dv[n, j] +:= 1; end if; end for; end for;
Im := RowSpace(BasisMatrix(K) * Dv);     // divisors (sums over the pole indices, i.e. Stab(d) * multiplicity)
printf "realisable divisors on S5 (coordinates Stab(d)*mult, order %o): dimension %o, basis\n%o\n", S5, Dimension(Im), BasisMatrix(Im);
// in multiplicities: divide coordinate d by Stab(d)
St := DiagonalMatrix(Rationals(), [ 1/Stab(d, D, N) : d in S5 ]);
printf "in multiplicities: \n%o\n", BasisMatrix(Im) * St;
// the sharp test: W = cusp forms vanishing at the empty indices; W_f = a_7 = -4 eigenspace of T(49)
empty := [ n : n in [1..B] | IsEmpty(cmd[n]) ];
KE := Kernel(Matrix(Rationals(), dim, #empty, [ A0[i, n] : n in empty, i in [1..dim] ]));
Wb := [ Vector(k) * A0 : k in Basis(KE) ];
printf "dim W_B = %o\n", #Wb;
// T(49) on W_B: a(49 n) + (-n/7) a(n) + 7 a(n/49), identified on admissible n <= B/49
p := 7; Bp := B div (p^2);
function Tp(g) return [ g[p^2*n] + KroneckerSymbol(-n, p) * g[n] + (n mod p^2 eq 0 select p * g[n div p^2] else 0) : n in [1..Bp] ]; end function;
Wmat := Matrix(Rationals(), [ [ g[n] : n in [1..Bp] ] : g in Wb ]);
Tmat := Matrix(Rationals(), [ Tp(g) : g in Wb ]);
ok, Tm := IsConsistent(Wmat, Tmat);   // Tm * Wmat = Tmat
printf "T(49) on W_B (consistent %o): char poly %o\n", ok, CharacteristicPolynomial(Tm);
Ef := Kernel(Tm + 4 * IdentityMatrix(Rationals(), #Wb));
Wf := [ &+[ v[i] * Wb[i] : i in [1..#Wb] ] : v in Basis(Ef) ];
printf "dim W_f = %o\n", #Wf;
// torsion test of a divisor on S5: exists c with the given divisor sums and c ⊥ W_f (c arbitrary elsewhere? no:
// the divisor outside S5 must be zero too)
function Torsion(mults, forms)
    rows := #forms + #cmds;
    A2 := ZeroMatrix(Rationals(), rows, B); rhs := ZeroMatrix(Rationals(), rows, 1);
    for i->g in forms do for n in [1..B] do A2[i, n] := g[n]; end for; end for;
    for j->d in cmds do
        for n in [1..B] do if d in cmd[n] then A2[#forms + j, n] +:= 1; end if; end for;
        k := Index(S5, d); rhs[#forms + j, 1] := k gt 0 select Stab(d, D, N) * mults[k] else 0;
    end for;
    return IsConsistent(Transpose(A2), Transpose(rhs));
end function;
cands := [ <"Delta = P-228 + P-760 - P-12 - P-40", [1, 1, -1, -1, 0]>,
           <"P-228 + P-12 - P-760 - P-40", [1, -1, 1, -1, 0]>,
           <"P-228 + P-40 - P-760 - P-12", [1, -1, -1, 1, 0]>,
           <"P-228 - P-760", [1, -1, 0, 0, 0]>, <"P-12 - P-40", [0, 0, 1, -1, 0]>,
           <"P-532 - P-228", [-1, 0, 0, 0, 1]>, <"P-532 - P-760", [0, -1, 0, 0, 1]>, <"P-532 - P-12", [0, 0, -1, 0, 1]>, <"P-532 - P-40", [0, 0, 0, -1, 1]>,
           <"2 P-532 - P-228 - P-760", [-1, -1, 0, 0, 2]>, <"2 P-532 - P-12 - P-40", [0, 0, -1, -1, 2]> ];
for t in cands do
    printf "%-40o torsion (W_f): %o | realisable (all of S): %o\n", t[1], Torsion(t[2], Wf), Torsion(t[2], [ [A0[i, n] : n in [1..B]] : i in [1..dim] ]);
end for;
quit;
