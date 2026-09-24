// The obstruction Eisenstein series is for the DUAL Weil representation, so its quadratic space is
// the NEGATIVE of the one the m>0 machinery uses: kappa*m = -r, not +r.  Note the code is already
// inconsistent about this -- its class number h/w is taken from Q(sqrt(-r)) while its Whittaker is
// taken from Q(sqrt(+r)).  This tests the consistent choice.
AttachSpec("ShimuraQuotients.spec");

btrue := AssociativeArray();
btrue[[15,2]] := [<2, 4>, <10, -4>, <30, 0>];
btrue[[6,5]]  := [<10, -6>, <15, -2>, <30, -6>];
btrue[[10,3]] := [<3, -2>, <12, -2>, <30, -6>];

for b in [[15,2],[6,5],[10,3]] do
    D := b[1]; N := b[2];
    Ld := ShimuraCurveLattice(D, N);
    detprimes := Set(PrimeDivisors(Determinant(ChangeRing(Ld`Q, Integers()))));
    printf "\n===== X0^%o(%o)  (kappa*m = -r) =====\n", D, N;
    for t in btrue[b] do
        r := Rationals()!t[1]; bt := Rationals()!t[2];
        km := -r;
        dd := Discriminant(Integers(QuadraticField(-(Numerator(r)*Denominator(r)))));
        chi := KroneckerCharacter(dd);
        K := QuadraticField(-(Numerator(r)*Denominator(r)));
        h := ClassNumber(K); wr := #TorsionSubgroup(UnitGroup(K));
        _, ch := IsSquare(Rationals()!(r/AbsoluteValue(dd)));
        Sc := Sort([p : p in detprimes join Set(PrimeDivisors(Numerator(r)))]);
        en := &*[Rationals() | 1 - Evaluate(chi, p)/p : p in Sc];
        ed := &*[Rationals() | 1 - 1/(Rationals()!p)^2 : p in Sc];
        gW := Rationals()!1; gB := Rationals()!1; locs := [];
        for p in Sc do
            if N mod p eq 0 then
                w := KYWhittaker54(km, p, 1);
                bp := w;
            else
                w := KYWhittaker53(km, p, D mod p eq 0);
                bp := KYbp(km, p, D mod p eq 0, 1/(Rationals()!p));
            end if;
            Append(~locs, <p, w>);
            gW *:= w; gB *:= bp;
        end for;
        V1 := ch*(Rationals()!h/wr)*(en/ed)*gW;
        V2 := ch*(Rationals()!h/wr)*gB;
        V3 := gW;
        V4 := ch*(Rationals()!h/wr)*gW;
        printf "  r=%-4o b_true=%-5o locs=%o\n", t[1], bt, locs;
        printf "        V1=%-12o V2=%-12o V3=%-12o V4=%-12o\n", V1, V2, V3, V4;
        printf "        ratios: %o\n",
               [ v eq 0 select "0" else Sprintf("%o", bt/v) : v in [V1,V2,V3,V4] ];
    end for;
end for;
quit;
