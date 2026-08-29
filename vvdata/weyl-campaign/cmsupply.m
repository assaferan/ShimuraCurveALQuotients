// CM-point SUPPLY triage, FAITHFUL version.
//
//   NORMALIZ_BIN=... magma -b DD:=58 NN:=3 cmsupply2.m
//
// The first attempt (cmsupply.m) measured RAW supply from RationalandQuadraticCMPoints and
// called almost everything FILTER_BOUND -- including 22_3 and 15_2, which both have working
// models.  It was measuring the wrong thing: production passes Keep := Set(d_divs)
// (EquationsCovers.m:225), so the hauptmodul anchors are exempt from the coprime filter, and
// the demand is not the MaxNum default 7 but num_vals = max(2g+5) over the target covers
// (EquationsCovers.m:237).
//
// This replicates the production chain exactly up to -- but not including -- the expensive
// value computation, reproducing the supply arithmetic of AbsoluteValuesAtCMPoints
// (SchoferFormula.m:1085-1118), whose shortfall raises "Could not find enough points".
// It needs the Borcherds forms (for d_divs), so it costs about what ppint.m costs.
//
// Output:
//   BASE D N demand <num_vals> genera <list>
//   POOL rat <n> quad <n> include <n>
//   TOPUP need <n> after_rat <n> after_fetch <n> after_quad <n>
//   CMVERD D N <OK|SHORT> margin <n>
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
D := 58; N := 3;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;

Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};

fs := BorcherdsForms(star, curves : Prec := 100);
d_divs := &cat[[T[1] : T in DivisorOfBorcherdsForm(f, star)] : f in [fs[-1], fs[-2]]];
Include := Set(d_divs);

// demand, exactly as EquationsCovers.m computes it (empty Targets => all covers)
genus_list := [curves[i]`g : i in star`CoveredBy];
num_vals := Maximum([2*g + 5 : g in genus_list]);
printf "BASE %o %o demand %o genera %o\n", D, N, num_vals, Sort(genus_list);

all_cm_pts := CandidateDiscriminants(star, curves : Keep := Include);
cm_pts_rat := all_cm_pts[1]; cm_pts_quad := all_cm_pts[2];
printf "POOL rat %o quad %o include %o\n", #cm_pts_rat, #cm_pts_quad, #Include;

// --- SchoferFormula.m:1085-1118, replicated ---
other_cm_rat  := [p : p in cm_pts_rat  | p[1] notin Include];
other_cm_quad := [p : p in cm_pts_quad | p[1] notin Include];
need := num_vals - #Include;
after_rat := 0; after_fetch := 0; after_quad := 0; short := false;
if #other_cm_rat ge need then
    after_rat := need; after_fetch := need; after_quad := need;
else
    after_rat := #other_cm_rat;
    need2 := need - #other_cm_rat;
    Excl := {p[1] : p in other_cm_rat} join Include;
    bd := 16;   // Maximum(include_bd*2, 16)
    new_rat, new_quad := RationalandQuadraticCMPoints(star : bd := bd, Exclude := Excl,
                                                      coprime_to_level := true, target := num_vals + 8);
    after_fetch := after_rat + #new_rat;
    need2 := need2 - #new_rat;
    if need2 gt 0 then
        avail_quad := #other_cm_quad + #new_quad;
        after_quad := after_fetch + Minimum(avail_quad, need2);
        if avail_quad lt need2 then short := true; end if;
    else
        after_quad := after_fetch;
    end if;
end if;
printf "TOPUP need %o after_rat %o after_fetch %o after_quad %o\n",
       need, after_rat, after_fetch, after_quad;
printf "CMVERD %o %o %o margin %o\n", D, N,
       short select "SHORT" else "OK", after_quad - need;
printf "DONE\n";
quit;
