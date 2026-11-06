// DualTest.m
//
// The dual test: a degree-zero CM divisor Delta on a quotient X of X_0^D(N) is torsion in its
// Jacobian as soon as some principal part with divisor Delta pairs to zero with a space of cusp
// forms of weight 3/2 that contains every component of the Heegner series of X
// (paper/hyperelliptic-15-38.tex, Section 2; Borcherds' modularity).  The test is sufficient only:
// a target it refuses is not shown to be of infinite order, and the precision B of the basis must
// reach every index the target needs (B >= |d|/4).
//
// Three variants, all at the scalar level M = 4 * oddpart(DN):
//   "sharp"  X = X_0^D(N)*: the space W of cusp forms vanishing at every index n where the Heegner
//            series Z(n) = sum_r (1/e_d) P_d, d = -4n/r^2, has no CM point;
//   "plus2"  2 || N, D odd: V = (U(4) + 2) W+, where W+ vanishes where the plus-type series
//            H(n) = sum_s P'_{-n/s^2} of X_0^D(N/2)* does (paper/hyperelliptic-21-58.tex); read to
//            precision B/4 against the dictionary of X;
//   "odd"    2 || N, D odd: Delta is pushed to X_0^D(N/2)* by the correspondence T,
//            T C(d) = (1 + chi_d(2)) C'(d) + (2 - chi_{d/4}(2)) C'(d/4), and T Delta is tested
//            against W+ in the plus-type dictionary (paper/hyperelliptic-39-34.tex).

declare verbose DualTest, 2;

DualTestData := recformat< D : RngIntElt, N : RngIntElt, Variant : MonStgElt, Prec : RngIntElt,
                           Space : SeqEnum, Cols : SeqEnum, Dict : Assoc, Discs : SeqEnum >;

intrinsic CarriesCMPoint(D::RngIntElt, N::RngIntElt, d::RngIntElt) -> BoolElt
{True if X_0^D(N) has CM points by the order of discriminant d (d < 0): every prime of D is
 non-split in Q(sqrt(d)) and prime to the conductor, and every prime of N not dividing the
 conductor is not inert.}
    d0 := FundamentalDiscriminant(d); _, f := IsSquare(d div d0);
    if exists{p : p in PrimeDivisors(D) | f mod p eq 0 or KroneckerSymbol(d0, p) eq 1} then return false; end if;
    if exists{p : p in PrimeDivisors(N) | f mod p ne 0 and KroneckerSymbol(d0, p) eq -1} then return false; end if;
    return true;
end intrinsic;

intrinsic CMPointWeight(D::RngIntElt, N::RngIntElt, d::RngIntElt) -> RngIntElt
{The order e_d of the stabiliser of a CM point of discriminant d on the star X_0^D(N)*: the unit
 group of the order, doubled when the point is fixed by an Atkin-Lehner involution w_m
 (d = -4m, or d = -m with m = 3 mod 4).  The Heegner series weights P_d by 1/e_d.}
    e := d eq -3 select 3 else (d eq -4 select 2 else 1);
    if exists{m : m in Divisors(D*N) | m gt 1 and GCD(m, (D*N) div m) eq 1 and (d eq -4*m or (m mod 4 eq 3 and d eq -m))} then
        e *:= 2;
    end if;
    return e;
end intrinsic;

// d = -4n/r^2 and d = -n/s^2 (the plus-type index set)
heeg_discs := func< n | [ -4*n div (r^2) : r in [1..Isqrt(4*n)] | (4*n) mod (r^2) eq 0 and ((-4*n div (r^2)) mod 4) in {0, 1} ] >;
plus_discs := func< n | [ -n div (s^2) : s in [1..Isqrt(n)] | n mod (s^2) eq 0 and ((-n div (s^2)) mod 4) in {0, 1} ] >;

function chi2(d)
    d0 := FundamentalDiscriminant(d); _, f := IsSquare(d div d0);
    return f mod 2 eq 0 select 0 else KroneckerSymbol(d0, 2);
end function;

intrinsic DualTestBasisDirectory() -> MonStgElt
{The directory of saved weight-3/2 bases: the environment variable DUALW_CACHE, else ~/dualw.}
    dir := GetEnv("DUALW_CACHE");
    if dir eq "" then dir := GetEnv("HOME") cat "/dualw"; end if;
    return dir;
end intrinsic;

intrinsic ReadWeightThreeHalvesBasis(M::RngIntElt, B::RngIntElt : Dir := DualTestBasisDirectory(), MD5 := "", Compute := false) -> Mtrx
{The cusp forms of weight 3/2 on Gamma_0(M) with the theta^3 multiplier, one per row, coefficients
 1..B, from the file Dir/basis_<M>_<B>.txt that tools/proof-*/basisS.m writes.  When MD5 is given
 the file must have that checksum.  With Compute, a missing file is computed and written (minutes
 at levels near 2000, days near 8000).}
    fn := Sprintf("%o/basis_%o_%o.txt", Dir, M, B);
    if not FileExists(fn) and Compute then
        vprintf DualTest, 1 : "computing the weight-3/2 cusp forms at level %o to precision %o\n", M, B;
        S := CuspidalSubspace(HalfIntegralWeightForms(Gamma0(M), 3/2));
        Bq := Basis(S, B + 2);
        A := Matrix(Rationals(), #Bq, B, [ Coefficient(g, n) : n in [1..B], g in Bq ]);
        Write(fn, Sprint(A, "Magma") : Overwrite := true);
    end if;
    require FileExists(fn) : Sprintf("no saved basis %o; compute it with tools/proof-4130-1/basisS.m", fn);
    if MD5 ne "" then
        sum := Split(Pipe(Sprintf("(md5sum %o 2>/dev/null || md5 -q %o) | cut -c1-32", fn, fn), ""))[1];
        require sum eq MD5 : Sprintf("%o has checksum %o, expected %o", fn, sum, MD5);
    end if;
    A0 := eval Read(fn);
    require Ncols(A0) ge B : Sprintf("%o has only %o coefficients", fn, Ncols(A0));
    return Submatrix(A0, 1, 1, Nrows(A0), B);
end intrinsic;

intrinsic DualTestSpace(D::RngIntElt, N::RngIntElt, A0::Mtrx : Variant := "sharp") -> Rec
{The test space of the dual test on X_0^D(N)* from the basis A0 of weight-3/2 cusp forms at level
 4 * oddpart(DN) (rows, coefficients 1..B), with the dictionary the targets are read in.  Variant
 is "sharp", "plus2" or "odd" (see the header of DualTest.m).}
    require Variant in {"sharp", "plus2", "odd"} : "Variant must be sharp, plus2 or odd";
    B := Ncols(A0); r := Nrows(A0);
    T := rec< DualTestData | D := D, N := N, Variant := Variant >;
    if Variant eq "sharp" then
        dict := AssociativeArray();
        for n in [1..B] do dict[n] := [ d : d in heeg_discs(n) | CarriesCMPoint(D, N, d) ]; end for;
        empty := [ n : n in [1..B] | IsEmpty(dict[n]) ];
        K := Kernel(Matrix(Rationals(), r, #empty, [ A0[i, n] : n in empty, i in [1..r] ]));
        T`Space := [ Eltseq(Vector(k) * A0) : k in Basis(K) ];
        T`Prec := B; T`Cols := [1..B]; T`Dict := dict;
    else
        require IsOdd(D) and N mod 2 eq 0 and (N div 2) mod 2 eq 1 : "plus2 and odd need D odd and 2 || N";
        Np := N div 2;
        Hzero := [ n : n in [1..B] | (n mod 4) in {1, 2} or forall{d : d in plus_discs(n) | not CarriesCMPoint(D, Np, d)} ];
        K := Kernel(Matrix(Rationals(), r, #Hzero, [ A0[i, n] : n in Hzero, i in [1..r] ]));
        Wp := [ Eltseq(Vector(k) * A0) : k in Basis(K) ];
        if Variant eq "plus2" then
            B4 := B div 4;
            T`Space := [ [ g[4*n] + 2*g[n] : n in [1..B4] ] : g in Wp ];
            dict := AssociativeArray();
            for n in [1..B4] do dict[n] := [ d : d in heeg_discs(n) | CarriesCMPoint(D, N, d) ]; end for;
            // every index is a column: a principal part may use indices with no CM point, which add
            // nothing to the divisor but change the pairing
            T`Prec := B4; T`Cols := [1..B4]; T`Dict := dict;
        else
            cols := [ n : n in [1..B] | n notin Hzero ];
            dict := AssociativeArray();
            for n in cols do dict[n] := [ d : d in plus_discs(n) | CarriesCMPoint(D, Np, d) ]; end for;
            T`Space := Wp; T`Prec := B; T`Cols := cols; T`Dict := dict;
        end if;
    end if;
    T`Discs := Sort(Setseq(&join([{Integers()|}] cat [ Set(T`Dict[n]) : n in T`Cols ])));
    vprintf DualTest, 1 : "DualTestSpace(%o, %o, %o): dimension %o at precision %o, %o CM discriminants\n",
        D, N, Variant, #T`Space, T`Prec, #T`Discs;
    return T;
end intrinsic;

intrinsic DualTestHecke(T::Rec, p::RngIntElt) -> BoolElt, AlgMatElt
{T(p^2) on the test space of a sharp dual test, a(p^2 n) + (-n/p) a(n) + p a(n/p^2), when the
 coefficients up to Prec/p^2 separate the space; false otherwise (the operator is then not
 determined by the coefficients available).}
    require T`Variant eq "sharp" : "the Hecke operator is read on the sharp space";
    Bp := T`Prec div p^2;
    if #T`Space eq 0 or Bp eq 0 then return false, _; end if;
    Tp := func< g | [ g[p^2*n] + KroneckerSymbol(-n, p) * g[n] + (n mod p^2 eq 0 select p * g[n div p^2] else 0) : n in [1..Bp] ] >;
    Wmat := Matrix(Rationals(), [ [ g[n] : n in [1..Bp] ] : g in T`Space ]);
    if Rank(Wmat) ne #T`Space then return false, _; end if;
    ok, Tm := IsConsistent(Wmat, Matrix(Rationals(), [ Tp(g) : g in T`Space ]));
    if not ok then return false, _; end if;
    return true, Tm;
end intrinsic;

intrinsic RestrictDualTestSpace(T::Rec, p::RngIntElt, a::RngIntElt) -> Rec
{The test space cut down to the eigenspace of T(p^2) with eigenvalue a (the part attached to the
 newforms with a_p = a), for a sharp test whose T(p^2) is determined.}
    ok, Tm := DualTestHecke(T, p);
    require ok : "T(p^2) is not determined on this space";
    E := Kernel(Tm - a * IdentityMatrix(Rationals(), #T`Space));
    S := T;
    S`Space := [ [ &+[ v[i] * T`Space[i][n] : i in [1..#T`Space] ] : n in [1..T`Prec] ] : v in Basis(E) ];
    return S;
end intrinsic;

intrinsic CertifiesTorsion(T::Rec, target::SeqEnum) -> BoolElt
{True if the dual test T certifies that the degree-zero CM divisor target = [<d, m>, ...] (m the
 multiplicity of the point P_d of X_0^D(N)*) is torsion: some principal part whose divisor is
 target pairs to zero with the test space.  False means only that this space does not certify it.}
    D := T`D; N := T`N;
    rhs_of := AssociativeArray();
    if T`Variant eq "odd" then
        Np := N div 2;
        for t in target do
            d := t[1]; alpha := CMPointWeight(D, N, d) * t[2];
            require CarriesCMPoint(D, N, d) : Sprintf("%o carries no CM point on X_0^%o(%o)", d, D, N);
            if CarriesCMPoint(D, Np, d) then
                rhs_of[d] := (IsDefined(rhs_of, d) select rhs_of[d] else 0) + alpha * (1 + chi2(d));
            end if;
            if d mod 4 eq 0 and ((d div 4) mod 4) in {0, 1} and CarriesCMPoint(D, Np, d div 4) then
                rhs_of[d div 4] := (IsDefined(rhs_of, d div 4) select rhs_of[d div 4] else 0) + alpha * (2 - chi2(d div 4));
            end if;
        end for;
    else
        for t in target do
            rhs_of[t[1]] := (IsDefined(rhs_of, t[1]) select rhs_of[t[1]] else 0) + CMPointWeight(D, N, t[1]) * t[2];
        end for;
    end if;
    require forall{d : d in Keys(rhs_of) | d in T`Discs} :
        Sprintf("a discriminant of %o is beyond the precision %o", target, T`Prec);
    cols := T`Cols; pos := AssociativeArray(); for j->n in cols do pos[n] := j; end for;
    sp := T`Space; nd := #T`Discs;
    A := ZeroMatrix(Rationals(), #sp + nd, #cols); rhs := ZeroMatrix(Rationals(), #sp + nd, 1);
    for i->g in sp do
        for j->n in cols do A[i, j] := g[n]; end for;
    end for;
    dpos := AssociativeArray(); for k->d in T`Discs do dpos[d] := k; end for;
    for j->n in cols do
        for d in T`Dict[n] do A[#sp + dpos[d], j] +:= 1; end for;
    end for;
    for d in Keys(rhs_of) do rhs[#sp + dpos[d], 1] := rhs_of[d]; end for;
    return IsConsistent(Transpose(A), Transpose(rhs));
end intrinsic;
