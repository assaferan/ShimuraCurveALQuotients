// Serre-Stark build, phase 2a: the theta-basis CATALOGUE on A2 = disc(L(2)).
//
// For each cyclic nondegenerate orthogonal summand class (t, u) of A2
// (|Z| = 2t, q(g) = u/4t), compute the complement C = Z^perp and the
// invariants of the Weil representation rho_C (weight-0 vectors), both for
// rho and its dual. A basis member is theta_{Z} tensor v with v invariant;
// its e_0(A2)-exponents are (u t) j^2, so honest ledger scales need u = 1
// (square class). Report per class: count, C-structure, invariant dims, and
// the y_0-component of each invariant (the e_0-pairing weight).
//
//   magma -b DD:=15 NN:=2 ss3.m
AttachSpec("ShimuraQuotients.spec");

D := 15; N := 2;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;

Ld := ShimuraCurveLattice(D, N);
G := 2 * ChangeRing(Ld`Q, Rationals());       // the 2-rescaled Gram
Gi := G^(-1);
den := LCM([Denominator(x) : x in Eltseq(Gi)]);
Ldual2 := RSpaceWithBasis(ChangeRing(den*Gi, Integers()));
L2 := RSpaceWithBasis(ScalarMatrix(3, den));
A, to_disc := Ldual2 / L2;
AA := AbelianGroup(Moduli(A));
toA := func< a | &+[ A | e[i]*A.i : i in [1..Ngens(A)] ] where e := Eltseq(AA!a) >;
vecof := func< eta | ChangeRing(toA(eta)@@to_disc, Rationals()) >;
qval := function(eta)
    v := vecof(eta); r := (v*G, v)/(2*den^2); return r - Floor(r);
end function;
bval := function(eta, mu)
    v := vecof(eta); w := vecof(mu); r := (v*G, w)/(den^2); return r - Floor(r);
end function;
printf "BASE %o %o  |A2| = %o  invariants %o\n", D, N, #AA, Moduli(A);

// ---- enumerate cyclic nondegenerate summands, grouped by (t, u-class) ----
seen := {}; classes := AssociativeArray();
for g in AA do
    if IsIdentity(g) then continue; end if;
    Z := sub< AA | g >;
    if Z in seen then continue; end if;
    Include(~seen, Z);
    c := #Z;
    if IsOdd(c) then continue; end if;   // D_2t-type only: |Z| = 2t
    hval := [ Integers()!(c*bval(AA.i, g)) mod c : i in [1..Ngens(AA)] ];
    Zc := AbelianGroup([c]);
    h := hom< AA -> Zc | [ hval[i]*Zc.1 : i in [1..Ngens(AA)] ] >;
    P := Kernel(h);
    if #(Z meet P) ne 1 then continue; end if;
    t := c div 2;
    // unit class: q(g') = u/4t over generators; normalize u to the smallest rep
    us := Sort([ Integers()!(4*t*qval(k*g)) : k in [1..c-1] | GCD(k, c) eq 1 ]);
    u := us[1];
    key := <t, u>;
    if not IsDefined(classes, key) then classes[key] := [* *]; end if;
    Append(~classes[key], <Z, P, g>);
end for;

// ---- Weil representation of a finite quadratic module (sub P of AA) ----
// rho(T) e_g = e(q(g)) e_g ;  rho(S) e_g = sig8 / sqrt|C| * sum_d e(-(g,d)) e_d
// with sig8 = Milgram: |C|^(-1/2) sum e(q) = e(sig/8). Dual rep: q -> -q.
weil_invariants := function(P : dual := false)
    els := [ x : x in P ];
    n := #els;
    lvl := LCM([ Denominator(qval(x)) : x in els ] cat [8]);
    K<z> := CyclotomicField(lvl);
    e := func< r | z^(Integers()!(lvl*(r - Floor(r)))) >;
    sgn := dual select -1 else 1;
    // rho(S)e_g = e(-sig/8)/sqrt(n) * sum_d e(-sgn*(g,d)) e_d, and by Milgram
    // e(-sig/8)/sqrt(n) = (sum_x e(-sgn*q(x)))/n exactly (|Milgram sum|^2 = n).
    gs2 := &+[ e(-sgn*qval(x)) : x in els ];
    // invariants lie in ker(T-1) = span of isotropic elements; impose S v = v
    // there: for v supported on iso, need (gs2/n)*sum_g e(-sgn*(g,d)) v_g = v_d
    // for ALL d (zero at non-isotropic d). n x #iso system, small kernel.
    iso := [ x : x in els | qval(x) eq 0 ];
    k := #iso;
    M := Matrix(K, n, k, [ [ (gs2/n)*e(-sgn*bval(g, d)) - (g eq d select 1 else 0)
                             : g in iso ] : d in els ]);
    I1 := Kernel(Transpose(M));   // solutions v (length k) with M*v = 0
    zero_idx := Index(iso, P!0);
    vecs := [ Eltseq(b) : b in Basis(I1) ];
    return Dimension(I1), [ v[zero_idx] : v in vecs ], iso, vecs;
end function;

printf "\nCATALOGUE (t, u): count, |C|, C-invts, dim inv(rho_C), dim inv(rho*_C), y0s\n";
ks := Sort([ k : k in Keys(classes) ]);
for key in ks do
    t := key[1]; u := key[2];
    reps := classes[key];
    Z, P, g := Explode(reps[1]);
    dim1, y01 := weil_invariants(P);
    dim2, y02 := weil_invariants(P : dual := true);
    printf "CLASS t=%o u=%o  count=%o  |C|=%o C=%o  dim=%o y0=%o  dimdual=%o y0dual=%o\n",
        t, u, #reps, #P, Invariants(P), dim1, y01, dim2, y02;
end for;
quit;
