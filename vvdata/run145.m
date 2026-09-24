// X0^14(5): the missing N > 1 validation of the m = 0 multiplier, in one run.
//
//   (1) GROUND TRUTH, from the Hauptmodul consistency relation alone (instant, no oracle).
//   (2) the SHIPPED closed form m0_multiplier, on every form of the base.
//   (3) the numeric ORACLE (1/2) c_eta(0), on the forms where the two can be compared.
//
// Sizing (measured by poles145.m, NOT guessed):
//   form  -2   -1    9    10   11   12   13   14   15
//   pole  91   16   70    35   70   70   16   35   70      max|eta exponent| 24-34
// The ground truth lives on rows -1 and -2 (cmtables.m: s = key -1, stilde = key -2).  Row -2 has
// pole order 91 and would force NumSamples = 448, Prec = 406 -- some 20 hours.  Row -1 has pole
// order 16 and needs only K = 128, Prec = 200, and testing EITHER row against its ground truth is a
// complete test of the formula at N = 5.  Form 13 has the same pole order and rides along for free
// (no ground truth for it -- it is a prediction).
//
//   magma vvdata/run145.m
AttachSpec("ShimuraQuotients.spec");

D := 14; N := 5;
KS := 128; PREC := 200;
if assigned KK then KS := StringToInteger(KK); end if;
if assigned PR then PREC := StringToInteger(PR); end if;
RUNFORMS := [-1, 13];

// =============== (1) ground truth from the consistency relation ==================================
// The uncorrected pipeline values (vvdata/gy_example36_crash.log); scale = scale_tilde = 1.
ds     := [ -4, -35, -280, -11, -91, -84, -51 ];
degs   := [  1,   1,    1,   1,   1,   1,   2 ];
s      := [* Infinity(), 5/8,   55/128,  0, 1, 1/4, 1/8 *];
stilde := [* Infinity(), 35/8, 585/128,  1, 0, 3/4, 9/4 *];
scale_tilde := stilde[Index(s,0)];
scale       := s[Index(stilde,0)];

firing := [];
for i->d in ds do
    Append(~firing, not IsEmpty(PrimeDivisors(N div GCD(N, FundamentalDiscriminant(d)))));
end for;
printf "X0^%o(%o)  firing discs: %o   non-firing: %o\n", D, N,
       [ds[i] : i in [1..#ds] | firing[i]], [ds[i] : i in [1..#ds] | not firing[i]];
error if not (firing[Index(s,0)] and firing[Index(stilde,0)]),
    "a normalising point is non-firing -- the ground-truth argument does not apply";

// u = N^-mult(f_-1), v = N^-mult(f_-2); one linear equation per NON-firing rational disc.
rows := [];
for i->d in ds do
    if firing[i] or degs[i] ne 1 or s[i] cmpeq Infinity() then continue; end if;
    Append(~rows, [Rationals() | s[i]/scale, stilde[i]/scale_tilde]);
    printf "  d = %-6o : (%o) u + (%o) v = 1\n", d, s[i]/scale, stilde[i]/scale_tilde;
end for;
error if #rows lt 2, "fewer than two non-firing rational points -- (u,v) is undetermined";
A := Matrix(Rationals(), #rows, 2, &cat rows);
ok, x := IsConsistent(Transpose(A), Vector(Rationals(), [1 : r in rows]));
error if not ok, "the non-firing equations are inconsistent: no multiplier pair explains them";
u := x[1]; v := x[2];
for r in rows do error if r[1]*u + r[2]*v ne 1, "over-determined system not satisfied"; end for;

function exact_log(w, N)
    e := Valuation(Numerator(w), N) - Valuation(Denominator(w), N);
    error if w ne (Rationals()!N)^e, Sprintf("%o is not a power of %o", w, N);
    return -e;
end function;
truth := AssociativeArray();
truth[-1] := exact_log(u, N);  truth[-2] := exact_log(v, N);
printf "  => u = %o, v = %o\n", u, v;
printf "  GROUND TRUTH: mult(form -1) = %o, mult(form -2) = %o\n\n", truth[-1], truth[-2];

// =============== (2) the shipped closed form =====================================================
t0 := Cputime();
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
ks := Sort([k : k in Keys(fs)]);
Ld := ShimuraCurveLattice(D, N);
M := IsOdd(D*N) select 4*D*N else 2*D*N;
printf "forms built (%o s): keys = %o, M = %o, |G| = %o\n\n", Cputime(t0), ks, M, #Ld`disc_grp;

printf "%-6o %-10o %-18o %o\n", "form", "pole", "m0_multiplier", "vs ground truth";
for k in ks do
    foo := qExpansionAtoo(fs[k], 1);
    mult := M0Multiplier(foo, qExpansionAt0(fs[k], 1), D, N, Ld);
    tag := IsDefined(truth, k)
           select Sprintf("%o  %o", truth[k], mult eq truth[k] select "MATCH" else "*** MISMATCH ***")
           else "";
    printf "%-6o %-10o %-18o %o\n", k, -Valuation(foo), mult, tag;
end for;

// =============== (3) the numeric oracle ==========================================================
printf "\noracle on forms %o at Prec = %o, NumSamples = %o ...\n", RUNFORMS, PREC, KS;
t1 := Cputime();
consts, isoelts, errs := VVConstantTerms([fs[k] : k in RUNFORMS], Ld, M
                                         : Prec := PREC, NumSamples := KS);
printf "done (%o s).  #isotropic = %o (2N-1 = %o)\n\n", Cputime(t1), #isoelts, 2*N-1;
printf "%-6o %-14o %-22o %-22o %o\n", "form", "max|PPerr|", "c_0(0)", "(1/2)c_eta(0)", "vs ground truth";
for i->k in RUNFORMS do
    mult := Real(consts[i][2])/2;
    tag := IsDefined(truth, k) select Sprintf("%o", truth[k]) else "(prediction)";
    printf "%-6o %-14o %-22o %-22o %o\n", k, ChangePrecision(errs[i], 6),
           ChangePrecision(consts[i][1], 10), ChangePrecision(mult, 12), tag;
    printf "       c_eta(0) at all nonzero isotropic: %o\n",
           [ChangePrecision(c, 10) : c in consts[i][2..#consts[i]]];
end for;
printf "\nGATE: discard any row whose max|PPerr| is not small -- see memory m0-multiplier-solved.\n";
printf "TOTAL %o s\n", Cputime(t0);
quit;
