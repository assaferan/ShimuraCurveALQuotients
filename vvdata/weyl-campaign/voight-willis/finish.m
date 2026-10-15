// From the numerical branch points to an equation over Q, with the twist fixed by counting points.
//
//   magma -b IN:=out/branch_35_1_1.5..m finish.m < /dev/null
//
// THE PROBLEM THE INVARIANTS SOLVE.  The expansions give the branch points of X/W -> X* in a
// coordinate on the star that the computation chose for itself, and nothing says that coordinate
// is defined over Q.  Writing the equation in it would therefore give irrational coefficients even
// though the curve is rational.  But the invariants of a binary form -- the Igusa-Clebsch ones in
// genus 2, the j-invariant in genus 1 -- do not depend on the coordinate at all, and for a curve
// defined over Q they are rational numbers.  So they can be read off numerically, recognised as
// rationals, and the curve rebuilt over Q from them (Mestre's construction, which Magma provides).
//
// WHAT IS LEFT AFTER THAT IS EXACTLY THE TWIST.  Invariants see the curve over Qbar, so the
// rebuilt curve is right up to a quadratic twist and no further.  The twist is settled by counting
// points against the trace formula: the two twists differ by the sign of a_p at every prime where
// the twisting number is not a square, so a handful of primes picks one out, and the choice is
// then confirmed at every good prime below the cap.
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
// `load` takes a string literal, so the file named on the command line is read and evaluated.
dat := eval (Read(IN) cat "\nreturn <D, N, W, genus, atinf, prec, branch>;");
D := dat[1]; N := dat[2]; W := dat[3]; genus := dat[4];
atinf := dat[5]; prec := dat[6]; branch := dat[7];
CCb := ComplexField(prec);

R4 := RealField(4);
PQ<x> := PolynomialRing(Rationals());
PC<T> := PolynomialRing(CCb);
CAP := StringToInteger(assigned PCAP select PCAP else "150");

printf "=== X_0^%o(%o) / %o, genus %o, %o branch points%o (precision %o) ===\n",
    D, N, W, genus, #branch + (atinf select 1 else 0),
    atinf select " including infinity" else "", prec;

// The binary form whose roots are the branch points.  A branch point at infinity just means the
// form has degree 2g+1 rather than 2g+2, which HyperellipticCurve and the invariants both accept.
f := &*[ T - z : z in branch ];
error if Degree(f) ne 2*genus + (atinf select 1 else 2),
    "the number of branch points does not match the genus";

// How many digits of the branch points are believed right.  It is set by the TRUNCATION of the
// power series, not by the working precision: the series is cut where rho^NT is 10^-T, so T is
// the number of digits the branch points carry however many digits the arithmetic is done in.
// The invariants are of degree up to ten in the coefficients of the sextic, so a few more digits
// go there.  Distinguish the two: the working precision must be comfortably ABOVE 80 or the
// invariants come out meaningless (at 40 they are wrong in every digit, which is cancellation in
// I10 and nothing to do with the input), while what a recognition may rely on is ACC.
// ACCD is the number of digits the INVARIANTS are believed right to, which is several fewer than
// the branch points carry: at a truncation of 1e-28 the absolute Igusa-Clebsch invariants came
// out right to 24 digits in a controlled test, so about four digits go into the degree-10
// arithmetic.  Set it from the truncation target of the run that produced the branch points,
// minus a margin, and never from the working precision.
ACC := StringToInteger(assigned ACCD select ACCD else Sprint(Max(8, prec div 6)));

// Recognise a real number as a rational, with the CLASSICAL CERTAINTY CRITERION rather than a
// stability test.  A rational p/q is the only one with denominator at most q within eps of r as
// soon as q^2 eps < 1/2, so that is what is required; anything weaker accepts noise.  It does,
// and did: at a truncation of 1e-28 the first two invariants at X_0^35(1) came back as
// 4242472886187715/291357 and -382940844641465/3748273, both of the right SIZE and both wrong,
// where the truth is 1498959544523423744/102942875 and -139237160005402624/102942875.  A
// stability test over two consecutive denominator bounds passed them.
//
// The criterion also says how many digits a base needs: an invariant with denominator q and
// numerator p needs about 2 log10(q) + log10(p/q) digits, so the seventeen-digit numerators and
// nine-digit denominators here want a truncation target near 1e-36.
// The error is not guessed from the truncation target but MEASURED: the invariant has to be
// real, so its imaginary part is a direct estimate of the noise in it, and taking twenty times
// that (floored by the stated accuracy) makes the criterion self-calibrating.  Hand-setting it
// cost two runs, once too loose and once too tight.
function Recognise(z, acc)
    r := Real(z);
    scale := Abs(r) gt 1 select Abs(r) else Parent(r) ! 1;
    eps := Max(10^(-acc) * scale, 20 * Abs(Imaginary(z)));
    for e in [2..acc] do
        q := BestApproximation(r, 10^e);
        if Abs(r - q) lt eps and Denominator(q)^2 * eps lt 1/2 then
            return q, Abs(r - q)/scale, "";
        end if;
    end for;
    return 0, 1.0, "no rational of small enough denominator fits within the accuracy";
end function;

reconstructed := false;
C := 0;
if genus eq 2 then
    I := IgusaClebschInvariants(f);
    error if Abs(I[1]) lt 10^(-prec div 3) * Abs(I[2])^(1/2),
        "I2 vanishes, so the normalisation used here does not apply";
    // the absolute invariants: unchanged by a coordinate change on the star AND by rescaling the
    // form, so they are the numbers that must come out rational
    js := [ I[1]^5/I[4], I[1]^3*I[2]/I[4], I[1]^2*I[3]/I[4] ];
    qs := []; errs := [];
    for z in js do
        q, e := Recognise(z, ACC);
        Append(~qs, q); Append(~errs, e);
    end for;
    printf "absolute Igusa-Clebsch invariants, recognised as rationals (relative errors %o):\n",
        [ R4!e : e in errs ];
    for i in [1..3] do printf "   j%o = %o\n", i, qs[i]; end for;
    error if exists{ e : e in errs | e gt 0.5 }, "an invariant was not recognised as a rational";
    S := [ Rationals() | 1, qs[2]/qs[1], qs[3]/qs[1], 1/qs[1] ];
    C := HyperellipticCurveFromIgusaClebsch(S);
    reconstructed := true;
elif genus eq 1 then
    // four branch points; the j-invariant of the double cover, from their cross-ratio
    b := branch cat (atinf select [] else []);
    error if #branch ne 4, "a genus-1 double cover needs four branch points";
    l := ((branch[1]-branch[3])*(branch[2]-branch[4]))/((branch[1]-branch[4])*(branch[2]-branch[3]));
    jv := 256*(l^2-l+1)^3/(l^2*(l-1)^2);
    q, e := Recognise(jv, ACC);
    printf "j-invariant, recognised as a rational (relative error %o): j = %o\n", R4!e, q;
    error if e gt 0.5, "the j-invariant was not recognised as a rational";
    C := EllipticCurveFromjInvariant(q);
    reconstructed := true;
else
    printf "genus %o: no reconstruction from invariants is implemented here\n", genus;
end if;
error if not reconstructed, "nothing to verify";

// Mestre's construction returns a model with enormous rational coefficients, which the reduction
// routines refuse until it is integral.  Clearing the denominators of f by the square of their
// common multiple is an isomorphism (y is scaled with them), after which the model can be reduced.
// Mestre's construction returns coefficients of several hundred digits, and such a model has
// spurious bad reduction nearly everywhere, so point counts cannot even be taken from it until it
// is minimised.  Clearing the denominators of f by the square of their common multiple is an
// isomorphism (y scales with them); then the minimal Weierstrass model brings the discriminant
// down to the true one.
function Tidy(Cq)
    if Type(Cq) eq CrvEll then return MinimalModel(Cq); end if;
    f, h := HyperellipticPolynomials(Cq);
    g := f + h^2/4;
    d := LCM([ Denominator(c) : c in Coefficients(g) ] cat [1]);
    Cr := HyperellipticCurve(d^2*g);
    try Cr := ReducedMinimalWeierstrassModel(Cr); catch e
        try Cr := ReducedModel(SimplifiedModel(Cr)); catch e2 ; end try;
    end try;
    return Cr;
end function;
C := Tidy(C);
printf "curve over Q from the invariants (correct up to a quadratic twist):\n   %o\n", C;

// Points over F_{p^d}, or -1 where the model does not reduce to a curve of the right genus
// there.  Reducing the polynomial rather than base-changing the curve keeps a degree drop or a
// singular fibre from being a crash.
function Counts(Cq, p, d)
    F := GF(p^d);
    if Type(Cq) eq CrvEll then
        // a quadratic twist of an integral model need not be integral, and the minimal model is
        // also the one whose discriminant says where the reduction is really bad
        Cm := MinimalModel(Cq);
        if Valuation(Integers() ! Discriminant(Cm), p) gt 0 then return -1; end if;
        // BaseChange of an elliptic curve over Q to a finite field refuses here even on an
        // integral minimal model, so the reduction is built from the a-invariants instead
        return #Points(EllipticCurve([ F | a : a in aInvariants(Cm) ]));
    end if;
    fq, hq := HyperellipticPolynomials(Cq);
    g := ChangeRing(fq + hq^2/4, F);
    if Degree(g) lt 5 or Discriminant(g) eq 0 then return -1; end if;
    return #Points(HyperellipticCurve(g));
end function;

// ---- the twist, from point counts against the trace formula ----
X := CreateShimuraQuot(D, N, Set(W));
X`g := genus;
good := [ p : p in PrimesUpTo(CAP) | p ne 2 and (D*N) mod p ne 0 ];
target := AssociativeArray();
for p in good do target[p] := [ ComputePointsViaTrace(X, p, d) : d in [1,2] ]; end for;

// THE TWISTING NUMBER, SOLVED FOR RATHER THAN SEARCHED.  Mestre's construction returns an
// arbitrary member of the isomorphism class over Qbar, and the number relating it to the curve
// over Q is not small: at X_0^85(1) it is 288862917341177298307577.  So it cannot be found by
// trying candidates.  It can be solved for.  Twisting a genus-2 curve by c replaces a_p by -a_p,
// so comparing the point count of the reconstructed curve with the trace formula at a prime p
// reads off the Legendre symbol (c|p) directly.  And c can only be divisible by primes where the
// reconstructed curve or the true one has bad reduction, which is the primes of the discriminant
// together with those of 2DN.  Writing c as a product of those with exponents 0 or 1, each prime
// p gives one linear equation over F_2 in those exponents and in the sign, and the system is
// heavily overdetermined.  The solution is then confirmed by the point counts it was not fitted
// to, over F_{p^2} as well as F_p.
dsc := Integers() ! (Numerator(Discriminant(C)) * Denominator(Discriminant(C)));
tf, un := TrialDivision(AbsoluteValue(dsc), 10^7);
sp := { Integers() | f[1] : f in tf } join Set(PrimeDivisors(2*D*N));
// whatever trial division leaves is attacked with a bounded elliptic-curve and sieve effort; a
// prime of the twisting number that escapes both makes the linear system below inconsistent,
// which is reported rather than guessed around
for u in un do
    try
        for q in Factorisation(u : ECMLimit := 2000, MPQSLimit := 50) do Include(~sp, q[1]); end for;
    catch e
        printf "  warning: a %o-digit factor of the discriminant did not factor\n", #Sprint(u);
    end try;
end for;
supp := Sort(SetToSequence(sp));
printf "the twisting number can only be supported on %o primes (largest %o)\n",
    #supp, #supp eq 0 select 0 else supp[#supp];

// one equation per prime where the count determines the symbol
rows := []; rhs := []; used := [];
for p in good do
    n1 := Counts(C, p, 1);
    if n1 lt 0 then continue; end if;
    t1 := target[p][1];
    if n1 eq t1 and 2*(p+1) - n1 eq t1 then continue; end if;   // a_p = 0: no information
    if n1 eq t1 then b := 0;
    elif 2*(p+1) - n1 eq t1 then b := 1;
    else
        printf "  p = %o: neither the curve nor its twist gives the trace formula count "
             * "(%o against %o and %o) -- the reconstruction is wrong\n", p, t1, n1, 2*(p+1)-n1;
        continue;
    end if;
    Append(~rows, [ GF(2) | (LegendreSymbol(q, p) eq 1 select 0 else 1) : q in supp ]
                  cat [ GF(2) | (p mod 4 eq 1 select 0 else 1) ]);   // the sign of c
    Append(~rhs, GF(2) ! b);
    Append(~used, p);
end for;
printf "%o primes give a usable symbol, %o unknowns\n", #rows, #supp + 1;
error if #rows lt #supp + 4, "too few primes to determine the twisting number; raise the prime cap";
Mx := Matrix(GF(2), #rows, #supp+1, &cat rows);
okv, vsol := IsConsistent(Transpose(Mx), Vector(GF(2), rhs));
error if not okv, "the symbols are inconsistent: the reconstructed curve is not a twist of the right one";
csgn := (vsol[#supp+1] eq GF(2)!1) select -1 else 1;
c := csgn * (&*([Integers() | 1] cat [ supp[i] : i in [1..#supp] | vsol[i] eq GF(2)!1 ]));
printf "twisting number solved: c = %o (kernel dimension %o)\n", c, #supp + 1 - Rank(Mx);

Cfin := Tidy(QuadraticTwist(C, c));
printf "\nMODEL X_0^%o(%o)/%o :\n   %o\n", D, N, W, Cfin;
nchk := 0; nskip := 0;
for p in good do
    for d in [1,2] do
        n1 := Counts(Cfin, p, d);
        if n1 lt 0 then nskip +:= 1; continue; end if;
        error if n1 ne target[p][d],
            Sprintf("point count disagrees at p = %o, degree %o: %o against %o", p, d, n1, target[p][d]);
        nchk +:= 1;
    end for;
end for;
printf "primes of bad reduction skipped: %o of %o\n", nskip, 2*#good;
printf "VERIFIED genus %o (formula %o) and %o point counts over F_p and F_p^2 against the trace formula\n",
    Genus(Cfin) , GenusShimuraCurveQuotient(D, N, Set(W)), nchk;
error if Genus(Cfin) ne GenusShimuraCurveQuotient(D, N, Set(W)), "genus mismatch";

// ---- if we already store a model of this quotient, the two must be isomorphic over Q ----
fn := Sprintf("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/data/models/models_%o_%o.m", D, N);
if FileExists(fn) then
    models := eval (Read(fn) cat "\nreturn models;");
    if IsDefined(models, W) and #models[W] ne 0 and Type(models[W][1][2]) ne MonStgElt then
        e := models[W][1];
        Cold := HyperellipticCurve(e[2] + (#e ge 3 select e[3]^2/4 else 0));
        printf "STORED MODEL comparison: isomorphic over Q: %o\n", IsIsomorphic(Cold, Cfin);
    end if;
end if;
if assigned PUBF then
    Cpub := (genus eq 1) select EllipticCurve(eval PUBF) else HyperellipticCurve(eval PUBF);
    printf "PUBLISHED equation comparison: isomorphic over Q: %o\n", IsIsomorphic(Cpub, Cfin);
end if;
printf "FINISH DONE\n";
quit;
