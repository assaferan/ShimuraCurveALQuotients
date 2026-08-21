// THE PAYOFF TEST.  Replace, at the ramified primes p | D only, the pipeline's lattice Whittaker by
// Proposition 5.3's anisotropic value L_p * b_p^{(2.17)}, and ask whether the assembled coefficient
// is proportional to the Eisenstein coefficients b_true solved from the vector-valued oracle.
// Everything else is left alone: p | N uses the code's factor (pinned to Prop 5.4 by
// tests/KudlaYangLocal.m) and the good primes use the code's factor (just shown to equal Prop 5.3).
AttachSpec("ShimuraQuotients.spec");

btrue := AssociativeArray();
btrue[[15,2]] := [<2, 4>, <10, -4>, <30, 0>];
btrue[[6,5]]  := [<10, -6>, <15, -2>, <30, -6>];
btrue[[10,3]] := [<3, -2>, <12, -2>, <30, -6>];

for b in [[15,2],[6,5],[10,3]] do
    D := b[1]; N := b[2];
    Ld := ShimuraCurveLattice(D, N);
    negQ := -ChangeRing(Ld`Q, Integers());
    detprimes := Set(PrimeDivisors(Determinant(ChangeRing(Ld`Q, Integers()))));
    Lfull := RSpaceWithBasis(IdentityMatrix(Integers(), 3));
    zero := Vector(Rationals(), [0,0,0]);
    printf "\n===== X0^%o(%o) =====\n", D, N;
    printf "%-5o %-7o %-13o %-13o %-13o %-13o\n", "r", "b_true",
           "A: code", "B: 5.3 at p|D", "C: B, no zeta", "D: bp only";
    for t in btrue[b] do
        r := Rationals()!t[1]; bt := Rationals()!t[2];
        d0 := -(Numerator(r)*Denominator(r));
        K := QuadraticField(d0); dd := Discriminant(Integers(K));
        chi := KroneckerCharacter(dd);
        h := ClassNumber(K); wr := #TorsionSubgroup(UnitGroup(K));
        _, ch := IsSquare(Rationals()!(r/AbsoluteValue(dd)));
        Sc := Sort([p : p in detprimes join Set(PrimeDivisors(Numerator(r)))]);
        en := &*[Rationals() | 1 - Evaluate(chi, p)/p : p in Sc];
        ed := &*[Rationals() | 1 - 1/(Rationals()!p)^2 : p in Sc];
        pre := ch * (Rationals()!h/wr) * (en/ed);
        gA := Rationals()!1; gB := Rationals()!1; gC := Rationals()!1; gD := Rationals()!1;
        for p in Sc do
            wcode := LocalWhittakerAtOne(r, p, zero, Lfull, negQ);
            ram := (D mod p eq 0);
            w53 := KYWhittaker53(r, p, ram);
            gA *:= wcode;
            gB *:= ram select w53 else wcode;
            gC *:= ram select w53*(1 - 1/(Rationals()!p)^2) else wcode;
            gD *:= KYbp(r, p, ram, 1/(Rationals()!p));
        end for;
        vals := [pre*gA, pre*gB, pre*gC, (Rationals()!h/wr)*ch*gD];
        printf "%-5o %-7o %-13o %-13o %-13o %-13o\n", t[1], bt, vals[1], vals[2], vals[3], vals[4];
        printf "      ratios: %o\n",
               [ v eq 0 select "P=0" else Sprintf("%o", bt/v) : v in vals ];
    end for;
end for;
quit;
