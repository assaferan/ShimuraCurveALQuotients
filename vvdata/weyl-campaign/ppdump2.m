// Dump each form's PRINCIPAL PART in exactly the (eta-bucket, r) indexing that m0_multiplier sums
// over, alongside the code's own per-index weight A^code.  With the oracle-measured multipliers this
// makes the true weights SOLVABLE rather than remembered: mult_f = sum_idx c_f(idx) * A_idx is a
// linear system in the A_idx, one equation per form.
//
//   magma DD:=15 NN:=2 ppdump.m
//
// Rebuilding the targets this way matters because a summarised A_m carries a normalisation that
// cannot be checked from the outside; here every quantity is emitted in the code's own convention, so
// A^solved and A^code are directly comparable.
//
// Output lines:
//   FORM  D N key mult_groundtruth
//   COEF  D N key blk idx r coeff        (blk = "oo" or "0"; idx is m at oo, j at cusp 0)
//   ACODE D N blk idx r A_code           (the weight m0_multiplier attaches to that index)
AttachSpec("ShimuraQuotients.spec");

D := 15; N := 2;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;

// Oracle-measured multipliers (memory m0-multiplier-solved; sum over isotropic eta divided by 4(N-1)).
truth := AssociativeArray();
truth[[15,2]] := [<-2,2>, <-1,4>, <9,0>, <10,0>, <11,4>, <12,2>, <13,4>, <14,-2>, <15,2>];
truth[[6,5]]  := [<-2,0>, <-1,0>, <9,6>, <10,3>, <11,12>, <12,9>, <13,3>, <14,0>, <15,3>];
truth[[10,3]] := [<-2,0>, <-1,0>, <9,6>, <10,3>, <11,3>, <12,3>, <13,3>, <14,0>, <15,6>];

Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
ks := Sort([k : k in Keys(fs)]);
Ld := ShimuraCurveLattice(D, N);
M := IsOdd(D*N) select 4*D*N else 2*D*N;
Q := Ld`Q; Qr := ChangeRing(Q, Rationals()); Qint := ChangeRing(Q, Integers());
negQ := -Qint; denom := Ld`denom; to_disc := Ld`to_disc;
detprimes := Set(PrimeDivisors(Determinant(Qint)));
Lfull := RSpaceWithBasis(IdentityMatrix(Integers(), 3));
printf "BASEINFO %o %o %o %o\n", D, N, M, Determinant(Qint);

// buckets, exactly as m0_multiplier forms them
mod_M_to_vecs := AssociativeArray([0..M-1]);
for j in [0..M-1] do mod_M_to_vecs[j] := []; end for;
i0 := 0;
for eta in Ld`disc_grp do
    if IsZero(eta) then i0 := eta; end if;
    v := ChangeRing(eta@@to_disc, Rationals());
    res := M*((v*Qr, v)/(2*denom^2));
    if not IsIntegral(res) then continue; end if;
    Append(~mod_M_to_vecs[Integers()!res mod M], eta);
end for;

// Emit, for one (eta, r), every ingredient separately -- the prefactor pieces and the per-prime local
// Whittaker VALUE and DERIVATIVE -- so that any candidate weight can be assembled on the python side
// without another Magma run.  The code's own weight is one such candidate, not a privileged one.
emit_contrib := procedure(blk, idx, eta, r, eta_id)
    w_eta := ChangeRing(eta@@to_disc, Rationals())/denom;
    D0 := -(Numerator(r)*Denominator(r));
    K := QuadraticField(D0);
    dd := Discriminant(Integers(K));
    chi := KroneckerCharacter(dd);
    h := ClassNumber(K); wr := #TorsionSubgroup(UnitGroup(K));
    ok, cond_half := IsSquare(Rationals()!(r/AbsoluteValue(dd)));
    error if not ok, "sqrt(r/|d|) irrational";
    Sc := Sort([p : p in detprimes join Set(PrimeDivisors(Numerator(r)))]);
    en := &*[Rationals() | 1 - Evaluate(chi, p)/p : p in Sc];
    ed := &*[Rationals() | 1 - 1/(Rationals()!p)^2 : p in Sc];
    printf "PREF %o %o %o %o %o %o %o %o %o %o %o\n", D, N, blk, idx, eta_id,
           r, cond_half, h, wr, en, ed;
    mu := Vector(Rationals(), Eltseq(w_eta));
    for p in Sc do
        Wq := LocalWhittakerPolynomial(r, p, mu, Lfull, negQ);
        ok2, Wr := CanChangeRing(Wq, Rationals());
        error if not ok2, "unscaled local polynomial is not rational";
        printf "LOCW %o %o %o %o %o %o %o %o\n", D, N, blk, idx, eta_id, p,
               Evaluate(Wr, 1), Evaluate(Derivative(Wr), 1);
    end for;
end procedure;

tr := AssociativeArray();
for t in truth[[D,N]] do tr[t[1]] := t[2]; end for;

oo_idxs := {};  zero_idxs := {};
for k in ks do
    printf "FORM %o %o %o %o\n", D, N, k, IsDefined(tr,k) select tr[k] else "?";
    foo := qExpansionAtoo(fs[k], 1);
    f0  := qExpansionAt0(fs[k], 1);
    for m in [1..-Valuation(foo)] do
        c := Coefficient(foo, -m);
        if c ne 0 then
            printf "COEF %o %o %o oo %o %o %o\n", D, N, k, m, m, c;
            Include(~oo_idxs, m);
        end if;
    end for;
    for j in [1..-Valuation(f0)] do
        c := Coefficient(f0, -j);
        if c ne 0 then
            printf "COEF %o %o %o 0 %o %o %o\n", D, N, k, j, (Rationals()!j)/M, c;
            Include(~zero_idxs, j);
        end if;
    end for;
end for;

for m in Sort(Setseq(oo_idxs)) do
    emit_contrib("oo", m, i0, Rationals()!m, 0);
end for;
for j in Sort(Setseq(zero_idxs)) do
    r := (Rationals()!j)/M;
    for i->eta in mod_M_to_vecs[j mod M] do
        emit_contrib("0", j, eta, r, i);
    end for;
end for;
printf "DONE\n";
quit;
