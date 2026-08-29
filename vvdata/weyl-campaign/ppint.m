// Integrality of the infinity-side principal parts, for EVERY Borcherds form on a base.
//
//   magma DD:=33 NN:=2 vvdata/weyl-campaign/ppint.m
//
// Schofer / Guo-Yang Thm B assume c_eta(m) in Z for m <= 0.  39_2 and 46_3 were both found
// to violate that (campaign 6c9c753, 9a34932).  The hypothesis under test is that
// NON-INTEGRALITY OF THE PRINCIPAL PART PREDICTS PIPELINE FAILURE -- i.e. the backlog
// blocker is form construction, not the m=0 term.  This dumps the evidence without running
// ValuesAtCMPoints, which is where the failing bases die.
//
// Schofer's hypothesis is c_eta(m) in Z for EVERY component, not just the one at infinity, so
// the 0-side principal part is checked too: a base that is oo-integral but 0-non-integral would
// otherwise be scored "integral" wrongly.  The oo verdict is reported separately because that is
// the block the earlier 39_2/46_3 measurements covered, and the two must stay comparable.
//
// Output:
//   BASE  D N #forms
//   PP    key blk m coeff             (blk = oo|0; one per nonzero principal-part coefficient)
//   VERD  key blk poleorder nterms denomlcm log10maxnum INTEGRAL|NONINTEGRAL
//   BASEVERD D N oo:INTEGRAL|NONINTEGRAL all:INTEGRAL|NONINTEGRAL
AttachSpec("ShimuraQuotients.spec");
D := 39; N := 2;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
ks := Sort([k : k in Keys(fs)]);
printf "BASE %o %o %o\n", D, N, #ks;
oo_ok := true; all_ok := true;
for k in ks do
    for blk in ["oo", "0"] do
        e := blk eq "oo" select qExpansionAtoo(fs[k], 1) else qExpansionAt0(fs[k], 1);
        if IsWeaklyZero(e) then continue; end if;
        v := Valuation(e);
        if v ge 0 then continue; end if;   // no principal part in this block
        cs := [];
        for n in [v..-1] do
            c := Coefficient(e, n);
            if c ne 0 then
                printf "PP %o %o %o %o\n", k, blk, -n, c;
                Append(~cs, Rationals()!c);
            end if;
        end for;
        dl := #cs eq 0 select 1 else LCM([Denominator(c) : c in cs]);
        mx := #cs eq 0 select 0 else Maximum([AbsoluteValue(Numerator(c)) : c in cs]);
        lg := mx eq 0 select 0 else Ilog(10, mx);
        ok := dl eq 1;
        if not ok then
            all_ok := false;
            if blk eq "oo" then oo_ok := false; end if;
        end if;
        printf "VERD %o %o %o %o %o %o %o\n", k, blk, -v, #cs, dl, lg,
               ok select "INTEGRAL" else "NONINTEGRAL";
    end for;
end for;
printf "BASEVERD %o %o oo:%o all:%o\n", D, N,
       oo_ok select "INTEGRAL" else "NONINTEGRAL",
       all_ok select "INTEGRAL" else "NONINTEGRAL";
printf "DONE\n";
quit;
