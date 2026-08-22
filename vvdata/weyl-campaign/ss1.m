// Serre-Stark census, step 1: which unary (cyclic, nondegenerate) orthogonal
// summands Z = <g> does A = L'/L admit, and what theta scales do they carry?
//
// Motivation: the ledger's weights w(N d k^2) must be 0-coset coefficients of a
// weight-1/2 theta for rho_A. Isotropic-quotient lifts alone cannot produce the
// observed scales (|A|/|B| is never a square at the needed |B|), so the thetas
// must sit on splittings A = Z + Z^perp, theta on Z tensor an invariant of Z^perp.
// For Z = disc(Z, t x^2) the e_0-exponents are t j^2 with t = c^2 q(g), c = |Z|.
//
//   magma -b DD:=15 NN:=2 ss1.m
AttachSpec("ShimuraQuotients.spec");

D := 15; N := 2;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;

Ld := ShimuraCurveLattice(D, N);
A := Ld`disc_grp; to_disc := Ld`to_disc; den := Ld`denom;
Qr := ChangeRing(Ld`Q, Rationals());
printf "BASE %o %o  |A| = %o  invariants %o\n", D, N, #A, Moduli(A);

AA := AbelianGroup(Moduli(A));
toA := func< a | &+[ A | e[i]*A.i : i in [1..Ngens(A)] ] where e := Eltseq(a) >;
vecof := func< eta | ChangeRing(toA(eta)@@to_disc, Rationals()) >;
qval := function(eta)  // Q(eta) mod 1, quadratic form
    v := vecof(eta);
    r := (v*Qr, v)/(2*den^2);
    return r - Floor(r);
end function;
bval := function(eta, mu)  // bilinear (eta, mu) mod 1
    v := vecof(eta); w := vecof(mu);
    r := (v*Qr, w)/(den^2);
    return r - Floor(r);
end function;

// enumerate cyclic subgroups via generators, dedupe by subgroup
seen := {};
summands := [* *];
for g in AA do
    if IsIdentity(g) then continue; end if;
    Z := sub< AA | g >;
    if Z in seen then continue; end if;
    Include(~seen, Z);
    c := #Z;
    // orthogonal complement: kernel of x |-> c * b(x, g) mod c
    hval := [ Integers()!(c*bval(AA.i, g)) mod c : i in [1..Ngens(AA)] ];
    Zc := AbelianGroup([c]);
    h := hom< AA -> Zc | [ hval[i]*Zc.1 : i in [1..Ngens(AA)] ] >;
    P := Kernel(h);
    // nondegenerate summand iff Z meet P = 1  (then A = Z + P, |P| = |A|/c)
    if #(Z meet P) ne 1 then continue; end if;
    // canonical q(g) over generators of Z
    qs := Sort([ qval(k*g) : k in [1..c-1] | GCD(k, c) eq 1 ]);
    scale := c^2 * qs[1];
    scales := Sort(Setseq({ c^2 * x : x in qs }));
    Append(~summands, < c, qs[1], scale, scales, Invariants(P), Z, P >);
end for;

printf "cyclic nondegenerate orthogonal summands: %o\n", #summands;
for s in summands do
    printf "  |Z| = %o  q(g) reps %o  scales c^2 q = %o  C = %o\n",
        s[1], s[4] eq [s[3]] select [s[2]] else s[4], s[4], s[5];
end for;

// which target scales (from the ledger: m-classes N*d) are available?
printf "LEDGER-RELEVANT: scales t with t = N*d, d squarefree, t <= %o:\n", D*N*4;
avail := {};
for s in summands do for t in s[4] do
    if Denominator(t) eq 1 then Include(~avail, Integers()!t); end if;
end for; end for;
printf "  integer scales available: %o\n", Sort(Setseq(avail));
quit;
