// The m = 0 multipliers at a COMPOSITE squarefree level, by support class (prop:composite of the
// standalone): M0MultipliersBySupport on X_0^6(35), M = 420, |L^v/L| = 88200, on eta quotients IN
// THE BORCHERDS INPUT SPACE -- weight 1/2 with the character of the pipeline's pool
// (lhs_integer_programming in BorcherdsForms.m: sum r = 1, sum d r = sum (M/d) r = 0 mod 24,
// prod d^r = 2 * square), found here as short vectors of the solution lattice; poles allowed at
// every cusp, so these are NOT Borcherds forms, only forms the coset sum is well defined for.
// ⚠ compmult.m .. compmult4.m used eta quotients OUTSIDE this space (weight 0, or weight 1/2 with the
// wrong character): the two-point and class-constancy checks then fail, correctly, since F_f is
// not even well defined on the cosets.
// What is checked: the routine runs at this size; the constant terms are constant on each support
// class {5}, {7}, {5,7}; whether the three classes carry DIFFERENT values.
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
SetVerbose("ShimuraQuotients", 2);
D := StringToInteger(D); N := StringToInteger(N);
M := IsOdd(D*N) select 4*D*N else 2*D*N;
Ld := ShimuraCurveLattice(D, N);
ds := Divisors(M); nd := #ds; ps := PrimeDivisors(M);
// unknowns: r (nd), then slack a, b (for the two mod-24 conditions) and c_p (one per prime)
nv := nd + 2 + #ps;
rows := [];
Append(~rows, [1 : d in ds] cat [0 : i in [1..2 + #ps]]);
Append(~rows, ds cat [-24, 0] cat [0 : p in ps]);
Append(~rows, [M div d : d in ds] cat [0, -24] cat [0 : p in ps]);
for i->p in ps do
    Append(~rows, [Valuation(d, p) : d in ds] cat [0, 0] cat [j eq i select -2 else 0 : j in [1..#ps]]);
end for;
A := Matrix(Integers(), rows);
rhs := Vector(Integers(), [1, 0, 0] cat [p eq 2 select 1 else 0 : p in ps]);
ok, x0 := IsConsistent(Transpose(A), rhs); assert ok;
K := KernelMatrix(Transpose(A));                 // integer kernel: rows span the solutions of A x = 0
Lk := LatticeWithBasis(K);
// shortest representatives of x0 + Lk, in the r-coordinates: reduce x0 against the kernel lattice
cl := Eltseq(ClosestVectors(Lk, -Vector(Rationals(), Eltseq(x0)) : Max := 1)[1]);
v0 := Eltseq(x0);
base := [Integers() | v0[i] + cl[i] : i in [1..nv]];
xs := [base];
// a few more solutions: add the shortest LLL-reduced kernel vectors (ShortVectors on this rank-24
// lattice enumerates far too many vectors and was killed for memory)
B := LLL(K);
for j in [1..3] do
    tv := Eltseq(B[j]);
    Append(~xs, [Integers() | base[i] + tv[i] : i in [1..nv]]);
end for;
rs := [x[1..nd] : x in xs];
for r in rs do
    assert &+r eq 1 and (&+[r[i]*ds[i] : i in [1..nd]]) mod 24 eq 0 and (&+[r[i]*(M div ds[i]) : i in [1..nd]]) mod 24 eq 0;
    assert IsSquare(&*[Rationals() | ds[i]^r[i] : i in [1..nd]] / 2);
    WriteStderr(Sprintf("r = %o   order at oo %o, at 0 %o\n", r, (&+[r[i]*ds[i] : i in [1..nd]]) div 24, (&+[r[i]*(M div ds[i]) : i in [1..nd]]) div 24));
end for;
R := EtaQuotientsRing(M, D*N);
fs := [EtaQuotient(R, r) : r in rs];
printf "X_0^%o(%o), M = %o, |L^v/L| = %o\n", D, N, M, #Ld`disc_grp;
t := Cputime();
arrs := M0MultipliersBySupport(fs, Ld, D, N);
printf "\n%o s\n", Cputime(t);
for i->a in arrs do
    printf "form %o: %o\n", i, [<p, a[p]> : p in Sort(SetToSequence(Keys(a)))];
end for;
