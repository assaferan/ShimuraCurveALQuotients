// AN EXACT MODEL OVER Q FOR A GENUS-TWO QUOTIENT -- ONE EVALUATION PER POINT.
//
//   magma -b DD:=35 NN:=1 WS:=1,5 RATDISCS:=-28,-35,-7 NC:=10 CUTA:=60 PRA:=60 \
//         PROG:=prog.txt tools/genus2model.m < /dev/null
//
// The method is written up in paper/nelson-method.tex, section "A genus-two quotient, and a
// model over Q".  The derivative is taken inside the lattice sum (NelsonDerivedRatio) rather than
// by Cauchy's formula on a circle, which replaces thirteen evaluations per point by one -- measured
// at 7.6 times faster -- and is the more accurate of the two, a circle of radius r amplifying the
// evaluation error by 1/r.
//
// THE CHAIN.  With omega_i = h_i dz the W-invariant differentials, x = h_2/h_1 is a hyperelliptic
// coordinate and y = (dx/dz)/h_1, so y^2 = F(x) is a linear system in the seven coefficients of F:
// seven sample points solve it and the rest test it.  The coordinate is then moved so that three
// rational CM points go to 0, 1 and infinity, in which coordinate every rational point of the
// curve has rational coordinate, so the sextic has rational coefficient ratios -- of small height,
// unlike the absolute invariants, whose numerators run to 10^17 and cannot be recognised at the
// precision available.  The remaining scalar is a quadratic twist, which tools/genus2twist.m reads off the
// point counts.
//
// A CM POINT AT AN ATKIN-LEHNER FIXED POINT IS HANDLED BY L'HOPITAL.  There h_1 and h_2 both
// vanish -- the invariant forms are pullbacks from the quotient and a pullback dies where the cover
// ramifies -- so x = h_2/h_1 is 0/0.  An Atkin-Lehner involution has order two, so both vanish to
// order exactly one and the limit is h_2'/h_1', which this route has already computed.  No circle
// is needed for these either.

AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
D := StringToInteger(DD); N := StringToInteger(NN);
W := { StringToInteger(x) : x in Split(WS, ",") };
rds := [ StringToInteger(x) : x in Split(RATDISCS, ",") ];
nc := StringToInteger(NC);
CUT := StringToInteger(CUTA); PR := StringToInteger(PRA);
R8 := RealField(8); C12 := ComplexField(12);
PQ<T> := PolynomialRing(Rationals());

procedure Note(s)
    Write(PROG, s);
end procedure;
Note(Sprintf("start: D = %o, N = %o, |W| = %o, cutoff %o, %o digits, %o sample points, %o CM points",
             D, N, #W, CUT, PR, nc, #rds));

function ALFixedDiscs(m)
    if m eq 1 then return {}; end if;
    if m eq 2 then return {-4, -8}; end if;
    if (m mod 4) eq 3 then return {-m, -4*m}; end if;
    return {-4*m};
end function;

B := NelsonBasis(D, N, W : Precision := PR, Cutoff := CUT);
error if B`g ne 2, "genus two only";
CC := B`CC; RR := RealField(Precision(CC)); A := B`A; O := B`ON;
bas := [ A ! x : x in Basis(O) ];
Lsp := Kernel(Transpose(Matrix(Integers(), [[ Integers() ! Trace(x) : x in bas ]])));
bl := [ &+[ b[i]*bas[i] : i in [1..4] ] : b in Basis(Lsp) ];
Qf := Matrix([[ Norm(x+y) - Norm(x) - Norm(y) : y in bl ] : x in bl ]);
Note("basis built");

function CMPointInH(d)
    found, v := FindLambda(Qf, -d, O, bl);
    if not found then return false, CC ! 0; end if;
    elt := &+[ v[i]/2*bl[i] : i in [1..#bl] ] + (((-d) mod 4 eq 3) select A ! (1/2) else A ! 0);
    th := A ! (2*elt - Trace(elt));
    if th^2 ne d then return false, CC ! 0; end if;
    M := ChangeRing(B`mrep(th), CC);
    a := M[1][1]; bb := M[1][2]; c := M[2][1];
    rt := Sqrt(CC ! (a^2 + bb*c));
    z := (a + rt)/c;
    if Imaginary(z) lt 0 then z := (a - rt)/c; end if;
    return true, z;
end function;

printf "=== X_0^%o(%o)/%o, genus 2, cutoff %o, %o digits\n", D, N, Sort(Setseq(W)), CUT, PR;

// --- the rational CM points, one evaluation each
cmx := []; cmkind := [];
for d in rds do
    ok, z := CMPointInH(d);
    error if not ok, Sprintf("the discriminant %o does not embed", d);
    h, hd, xv, _ := NelsonDerivedRatio(B, B`z1, z : Cutoff := CUT);
    br := exists{ m : m in W | d in ALFixedDiscs(m) };
    if br then
        // both forms vanish here; the coordinate is the ratio of their derivatives
        xv := hd[2]/hd[1];
    end if;
    Append(~cmx, xv); Append(~cmkind, br);
    Note(Sprintf("CM d = %-7o %-10o |h1| = %-14o x = %o", d,
         br select "fixed pt" else "ordinary", R8 ! Abs(h[1]), C12 ! xv));
end for;
printf "\n%-8o %-10o %-32o %o\n", "d", "kind", "x", "Im(x)/Re(x)";
for i in [1..#rds] do
    printf "%-8o %-10o %-32o %o\n", rds[i], cmkind[i] select "fixed pt" else "ordinary",
           C12 ! cmx[i], RealField(17) ! (Imaginary(cmx[i])/Real(cmx[i]));
end for;
rsp := [ Imaginary(t)/Real(t) : t in cmx ];
printf "spread of Im(x)/Re(x): %o   (the rational points lie on one line in this normalisation)\n",
       R8 ! (Maximum(rsp) - Minimum(rsp));

// --- the sextic, one evaluation per sample point
pts := NelsonSamplePoints(B, nc);
xs := [ CC | ]; ys := [ CC | ];
for j in [1..nc] do
    h, hd, xv, dxv := NelsonDerivedRatio(B, B`z1, pts[j] : Cutoff := CUT);
    Append(~xs, xv); Append(~ys, dxv/h[1]);
    Note(Sprintf("sample %o/%o done", j, nc));
end for;
nu := 7;
M := Matrix(CC, nu, nu, [ xs[i]^(k-1) : k in [1..nu], i in [1..nu] ]);
error if Abs(Determinant(M)) lt 10^(-PR div 2), "the seven sample points do not determine F";
cs := Eltseq(M^-1 * Matrix(CC, nu, 1, [ ys[i]^2 : i in [1..nu] ]));
PCx<xx> := PolynomialRing(CC);
F := &+[ cs[k]*xx^(k-1) : k in [1..nu] ];
worst := RR ! 0;
for i in [nu+1..nc] do
    worst := Max(worst, Abs(Evaluate(F, xs[i]) - ys[i]^2)/Abs(ys[i]^2));
end for;
printf "\nworst prediction at the %o sample points held back: %o\n", nc - nu, R8 ! worst;
error if nc le nu, "use more than seven sample points, so that at least one is held back";

// --- the coordinate sending the first three CM points to 0, 1, infinity
x1 := cmx[1]; x2 := cmx[2]; x3 := cmx[3];
Ac := x2 - x3; Bc := -x1*(x2 - x3); Cc := x2 - x1; Dc := -x3*(x2 - x1);
PCu<uu> := PolynomialRing(CC);
G := &+[ cs[k+1] * (-Dc*uu + Bc)^k * (Cc*uu - Ac)^(6-k) : k in [0..6] ];
gc := [ Coefficient(G, k) : k in [0..6] ];
_, imax := Maximum([ Abs(t) : t in gc ]);
gn := [ t/gc[imax] : t in gc ];
// The threshold is the accuracy the fit actually achieved, measured at the points held back, not
// the working precision: a coefficient that is genuinely zero shows up at that size, while a wrong
// recognition is off by order one.
tol := Maximum(100*worst, RR ! 10^(-PR + 5));
printf "\nthe sextic with coefficient %o scaled to 1; tolerance %o from the measured accuracy\n",
       imax-1, R8 ! tol;
Qb := Maximum(10, Floor(Sqrt(1/tol)/100));
printf "denominator bound %o, so a chance fit would only reach %o\n", Qb, R8 ! (RR ! 1/Qb^2);
printf "%-5o %-34o %-26o %-14o %o\n", "k", "coefficient", "nearest rational", "error", "better than chance by";
rat := true; qs := [];
for k in [1..7] do
    z := gn[k];
    q := BestApproximation(Real(z), Qb);
    err := Max(Abs(Imaginary(z)), Abs(RR ! q - Real(z)))/(1 + Abs(z));
    printf "%-5o %-34o %-26o %-14o %o\n", k-1, C12 ! z, q, R8 ! err,
           (err gt 0) select Sprint(R8 ! ((RR ! 1/Qb^2)/err)) else "exact";
    if err gt tol then rat := false; end if;
    Append(~qs, q);
end for;
if not rat then printf "\n=> not rational at this accuracy; raise the cutoff\n"; quit; end if;

den := LCM([ Denominator(q) : q in qs ]);
GQ := &+[ (Integers() ! (den*qs[k+1])) * T^k : k in [0..6] ];
printf "\nsextic over Q, up to a quadratic twist: %o\n", GQ;
Note(Sprintf("sextic: %o", GQ));
jj := 0; ok := true;
try jj := G2Invariants(HyperellipticCurve(GQ)); catch e ok := false; end try;
if not ok then printf "=> singular\n"; quit; end if;
printf "absolute invariants: %o\n", jj;

// --- the CM points that fixed nothing are a test: they too must be rational here
printf "\nthe other rational CM points, which fixed nothing:\n";
allrat := true;
for i in [4..#rds] do
    u := (Ac*cmx[i] + Bc)/(Cc*cmx[i] + Dc);
    q := BestApproximation(Real(u), Qb);
    err := Max(Abs(Imaginary(u)), Abs(RR ! q - Real(u)))/(1 + Abs(u));
    printf "   d = %-7o u = %-30o nearest rational %-16o error %o\n", rds[i], C12 ! u, q, R8 ! err;
    if err gt tol then allrat := false; end if;
end for;
printf "all of them rational: %o\n", allrat;

// --- the twist, which completes the model
//
// WHY A TWIST IS LEFT AT ALL, since the curve was read off its own differentials.  Nelson's formula
// returns g_i = y_1 y_2 conj(f_i(z_1)) f_i(z_2), so each form carries an unknown constant of the
// base point, h_i = K_i f_i.  The coordinate x = h_2/h_1 is only rescaled by that, and the Mobius
// step absorbs it; but y = (dx/dz)/h_1 is multiplied by mu = K_2/K_1^2, so the fitted sextic is
// mu^2 times the honest one.  Dividing by the largest coefficient above then discards mu^2
// outright.  So the scalar is LOST, not absent: |mu|^2 could be recovered by evaluating at
// z_2 = z_1, which returns y_1^2 |f_i(z_1)|^2, and the phase is pinned by the coefficients having
// to be real -- leaving only a SIGN, the twist by -1, which one prime settles.  The point counts
// below do that job, and are worth keeping anyway as an independent check of the isogeny class.
//
// ⚠ The counts fix the twist and confirm the isogeny class; they could NOT select the curve, since
// every principal polarization on the surface shares an L-function.  The curve comes from the
// differentials, not from here.
supp := Set(PrimeDivisors(D*N))
        join Set(PrimeDivisors(Integers() ! (den^(2*Degree(GQ)) * Discriminant(GQ))));
for c in Coefficients(GQ) do
    if c ne 0 then supp join:= Set(PrimeDivisors(Numerator(c))); end if;
    supp join:= Set(PrimeDivisors(Denominator(c)));
end for;
supp := Sort(SetToSequence(supp));
function CountAt(poly, p)
    dd := LCM([ Denominator(c) : c in Coefficients(poly) ] cat [1]);
    if dd mod p eq 0 then return -1; end if;
    gg := ChangeRing(PQ ! (dd^2 * poly), GF(p));
    if Degree(gg) lt Degree(PQ ! (dd^2*poly)) - 1 or Degree(gg) lt 3 or Discriminant(gg) eq 0 then
        return -1;
    end if;
    return #Points(HyperellipticCurve(gg));
end function;
cnt := AssociativeArray();
for p in PrimesInInterval(3, 300) do
    if (D*N) mod p ne 0 then
        c := CountAt(GQ, p);
        if c ge 0 then cnt[p] := c; end if;
    end if;
end for;
gp := Sort([ p : p in Keys(cnt) ]);
MS := ModularSymbols(D*N, 2);
hit := 0; traces := AssociativeArray();
for Fm in NewformDecomposition(NewSubspace(CuspidalSubspace(MS))) do
    if Dimension(Fm) ne 4 then continue; end if;
    tr := AssociativeArray(); good := true;
    for p in gp do
        // Trace(HeckeOperator) is TWICE Tr(a_p): T_p acts on the four-dimensional homology as a_p
        t := (Integers() ! Trace(HeckeOperator(Fm, p))) div 2;
        tr[p] := t;
        // the orbit is matched on |Tr a_p|, which is free of the Atkin-Lehner sign convention
        if Abs(p + 1 - cnt[p]) ne Abs(t) then good := false; break; end if;
    end for;
    if good then hit := Fm; traces := tr; break; end if;
end for;
if hit cmpeq 0 then
    printf "\n=> no dimension-two newform orbit matches |Tr a_p| at all %o good primes, so this\n", #gp;
    printf "   curve is not in the isogeny class of any of them; the twist is not determined\n";
else
    printf "\nmatched a newform orbit on |Tr a_p| at all %o good primes\n", #gp;
    rows := []; rhs := []; used := [];
    for p in gp do
        t := traces[p];
        if t eq 0 then continue; end if;       // carries no information about the character
        b := (p + 1 - cnt[p] eq t) select 0 else 1;
        Append(~rows, [ GF(2) | (LegendreSymbol(q, p) eq 1 select 0 else 1) : q in supp ]
                      cat [ GF(2) | (p mod 4 eq 1 select 0 else 1) ]);
        Append(~rhs, GF(2) ! b); Append(~used, p);
    end for;
    if #used lt #supp + 4 then
        printf "=> too few informative primes to determine the twist\n";
    else
        Mt := Matrix(GF(2), #rows, #supp+1, &cat rows);
        ok2, vv := IsConsistent(Transpose(Mt), Vector(GF(2), rhs));
        if not ok2 then
            printf "=> the Legendre symbols are inconsistent: no single twist fits the counts\n";
        elif Rank(Mt) ne #supp+1 then
            printf "=> the twist is not determined uniquely; widen the prime range\n";
        else
            sgn := (vv[#supp+1] eq GF(2)!1) select -1 else 1;
            tw := sgn * (&*([Integers() | 1] cat [ supp[i] : i in [1..#supp] | vv[i] eq GF(2)!1 ]));
            printf "twisting constant: %o\n", tw;
            Cm := ReducedMinimalWeierstrassModel(IntegralModel(HyperellipticCurve(tw*GQ)));
            printf "\nTHE MODEL: %o\n", Cm;
            Note(Sprintf("model: %o", Cm));
            bad := 0;
            for p in gp do
                if CountAt(tw*GQ, p) ne p + 1 - traces[p] then bad +:= 1; end if;
            end for;
            printf "counts against the newform: %o mismatches over %o good primes\n", bad, #gp;
        end if;
    end if;
end if;

fn := Sprintf("data/models/models_%o_%o.m", D, N);
if FileExists(fn) then
    models := eval (Read(fn) cat "\nreturn models;");
    key := Sort([ Integers() | w : w in W ]);
    if IsDefined(models, key) and #models[key] gt 0 and models[key][1][1] eq 2 then
        e := models[key][1];
        want := G2Invariants(HyperellipticCurve(PQ ! e[2], PQ ! e[3]));
        printf "\nstored model's invariants: %o\n", want;
        printf "%o\n", (jj eq want) select "=> EXACT MATCH with the independently known model"
                       else "=> the invariants differ from the stored model";
    end if;
end if;
quit;
