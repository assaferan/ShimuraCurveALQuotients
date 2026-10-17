// THE GENUS-4 CRITERION, TESTED WITHOUT A MODEL OF THE STAR CURVE.
//
//   magma -b DISC:=2002 LEV:=1 MA:=11 MB:=91 PR:=60 CUT:=40 NP:=8 RHO:=100 g4test.m
//
// X = X_0^D(N)/W of genus 4 sits over the genus-2 star curve C by an Atkin-Lehner involution,
// branched at two points, and if X is hyperelliptic those two must have the SAME IMAGE under
// the canonical map of C (they are then exchanged by its hyperelliptic involution, i.e.
// P + Q ~ K_C).  The canonical map of C is the ratio of its two holomorphic differentials,
// which are the two W*-invariant weight-2 newforms of level D*N, so the whole test is a
// comparison of two numbers -- far more robust than a rank determination.
//
// TWO THINGS THAT HAVE TO BE HANDLED.
//
// (a) THE RATIO IS 0/0 AT THE BRANCH POINTS.  A w_m-invariant differential vanishes at every
//     fixed point of w_m: the involution acts by -1 on the cotangent space there, so an
//     invariant differential must vanish.  Concretely, for mu of trace 0 and reduced norm m the
//     fixed point is z0 = (a + i sqrt(m))/c and (c z0 + d)^2 = -m, so the automorphy relation
//     m (c z0 + d)^-2 f(z0) = eps f(z0) reads -f(z0) = eps f(z0), forcing f(z0) = 0 when
//     eps = +1.  Both differentials vanish there, so the ratio is analytic but indeterminate,
//     and its value is recovered as the MEAN OVER A SMALL CIRCLE (mean value property), with an
//     error of order rho^NP.  That the values at the centre really do collapse against the
//     circle is printed below and is a direct check that z0 is a genuine Atkin-Lehner fixed
//     point rather than an arbitrary point.
//
// (b) THE OVERALL CONSTANT CANCELS.  Theorem 3.1 returns Im(z1) Im(z) conj(f_i(z1)) f_i(z), so
//     the ratio of the two values carries an unknown factor conj(f_1(z1))/conj(f_2(z1)).  It is
//     the SAME at both branch points, so it cancels in the comparison, which is why no
//     normalisation and no recognition of the individual values is needed.
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/nelson/nelson3.m";

D := StringToInteger(DISC);
N := (assigned LEV) select StringToInteger(LEV) else 1;
ms := [ StringToInteger(MA), StringToInteger(MB) ];
ds := [ StringToInteger(DA), StringToInteger(DB) ];
PREC := StringToInteger(PR); cutv := StringToInteger(CUT);
np := StringToInteger(NP); frac := StringToInteger(RHO)/10000;
NREP := StringToInteger(assigned NR select NR else "1");
CC := ComplexField(PREC); RR := RealField(PREC); Ic := CC.1; R8 := RealField(10);

k := 2; DN := D*N;
full := { Integers() | d : d in Divisors(DN) | GCD(d, DN div d) eq 1 };
gC := GenusShimuraCurveQuotient(D, N, full);
printf "=== star curve of X_0^%o(%o): genus %o; testing the branch points w_%o and w_%o\n",
       D, N, gC, ms[1], ms[2];
require gC eq 2 : "this test is for a genus-2 star curve";

// ---- the two differentials of the star curve ----
nmax := Max(60, Ceiling(DN*(cutv/(4*Pi(RR)))^2) + 20);
NFs := Newforms(CuspForms(Gamma0(DN), k));
sel := [];
for orb in NFs do
    for f in orb do
        if Level(Parent(f)) ne DN then continue; end if;
        KF := BaseRing(Parent(f)); qq := qExpansion(f, nmax+1);
        if Degree(KF) eq 1 then
            cs := [ CC ! (Rationals() ! Coefficient(qq, n)) : n in [1..nmax] ];
        else
            cj := Conjugates(KF.1 : Precision := PREC+10); idx := Index(orb, f);
            cs := [ CC | Evaluate(Polynomial(Rationals(), Eltseq(Coefficient(qq,n))), cj[idx])
                    : n in [1..nmax] ];
        end if;
        inv := true;
        for m in full do
            if m eq 1 then continue; end if;
            e := CC ! 1;
            for p in PrimeDivisors(m) do
                e *:= (D mod p eq 0) select cs[p] else -cs[p];
            end for;
            if Abs(e - 1) gt 0.01 then inv := false; end if;
        end for;
        if inv then Append(~sel, cs); end if;
    end for;
end for;
printf "W*-invariant newforms of level exactly %o: %o (genus of the star says %o)\n",
       DN, #sel, gC;
require #sel eq gC : "the W*-invariant newforms of level exactly D*N do not span";

Wtab := WkTable(k, cutv + 2, RR);
dat := NelsonSetupT(D, N, k, PREC, cutv + 2, Wtab);
A := dat`A; O := dat`O; mrep := dat`mrep;

// ---- Atkin-Lehner fixed points, from OPTIMAL EMBEDDINGS ----
// A box or short-vector search is the wrong shape for an order and silently misses elements;
// Embed returns an optimal embedding of the quadratic order, so the fixed point really has the
// discriminant asked for and not that of some larger order containing it.  For d = 1 mod 4 the
// element sqrt(d) has trace 0 and reduced norm |d| = m, which is the Atkin-Lehner element of
// w_m, and its fixed point in the upper half plane is (a + i sqrt(m))/c.
pts := [];
for j in [1..#ms] do
    m := ms[j]; d := ds[j];
    K := QuadraticField(d);
    R := MaximalOrder(K);
    require Discriminant(R) eq d : "the quadratic order does not have the branch discriminant";
    _, emb := Embed(R, O);
    mu := A ! emb(R ! (K.1));
    require Trace(mu) eq 0 and Norm(mu) eq m :
        "the embedded square root is not the Atkin-Lehner element of w_m";
    require forall{ b : b in ZBasis(O) | (-mu*(A!b)*mu/m) in O } :
        "the element does not normalise the order";
    M := ChangeRing(mrep(mu), RR);
    a := M[1][1]; c := M[2][1];
    z := (a + Sqrt(RR!m)*Ic)/c;
    if Imaginary(z) lt 0 then z := (a - Sqrt(RR!m)*Ic)/c; end if;
    printf "w_%o (discriminant %o): fixed point z0 = %o, Im = %o\n", m, d,
        ComplexField(14)!z, RealField(10)!Imaginary(z);
    Append(~pts, <m, d, z>);
end for;

// ---- evaluate ----
z1 := CC ! (0.1103 + 1.0411*Ic);
om := Exp(2*Pi(CC)*Ic/np);
results := [];
t0 := Cputime();
for pp in pts do
    m := pp[1]; z0 := pp[3];
    // the circle must be small compared with the height, or it reaches the real axis
    rho := (RR ! frac) * Imaginary(z0);
    Sc, _ := NelsonEvalMulti(dat, z1, z0, sel, RR!cutv, CC);
    rs := [ CC | ]; mags := [ RR | ];
    for t in [0..np-1] do
        zz := z0 + rho*om^t;
        S, _ := NelsonEvalMulti(dat, z1, zz, sel, RR!cutv, CC);
        Append(~rs, S[1]/S[2]);
        Append(~mags, Abs(S[2])/Imaginary(zz));
    end for;
    circmag := &+mags/np;
    cenmag := Abs(Sc[2])/Imaginary(z0);
    rval := (&+rs)/np;
    rhalf := (&+[ rs[2*t+1] : t in [0..(np div 2)-1] ])/(np div 2);
    printf "\nw_%o at z0 = %o  (circle radius %o)\n", m, ComplexField(14)!z0, RealField(8)!rho;
    printf "   |f2| at the centre %o against %o on the circle: ratio %o\n",
        R8!cenmag, R8!circmag, R8!(cenmag/circmag);
    printf "   (a W*-invariant differential MUST vanish at an Atkin-Lehner fixed point, so a "
         * "small ratio here is the check that z0 really is one)\n";
    printf "   canonical-map value (circle mean of f1/f2) = %o\n", ComplexField(24)!rval;
    printf "   same from half the circle points: %o   (relative difference %o)\n",
        ComplexField(14)!rhalf, R8!(Abs(rval-rhalf)/Abs(rval));
    Append(~results, <m, rval>);
end for;
printf "\ntotal evaluation time %o s\n", R8!Cputime(t0);

printf "\n==== THE TEST ====\n";
for i in [1..#results] do
    printf "  w_%-4o  canonical image %o\n", results[i][1], ComplexField(24)!results[i][2];
end for;
if #results ge 2 then
    v1 := results[1][2]; v2 := results[#results][2];
    d := Abs(v1-v2)/Maximum(Abs(v1), Abs(v2));
    printf "\n  relative difference between the two branch points: %o\n", R8!d;
    printf "  => %o\n", d lt 10^(-10)
        select "THE TWO BRANCH POINTS HAVE THE SAME IMAGE: P + Q ~ K_C, consistent with X hyperelliptic"
        else "THE TWO IMAGES DIFFER: P + Q is NOT the canonical class.  X is hyperelliptic only if BOTH points are Weierstrass points of the star curve, which is the remaining case to exclude";
end if;
printf "G4TEST DONE\n";
quit;
