// THE QUADRATIC TWIST OF A GENUS-TWO MODEL, FROM POINT COUNTS AGAINST THE NEWFORM.
//
//   magma -b DD:=35 NN:=1 SEXTIC:="-49/100*u^5 + u^4 - 277/350*u^3 + 69/245*u^2 - 1/4900*u" \
//         PRIMECAP:=200 tools/genus2twist.m < /dev/null
//
// The route of tools/genus2model.m normalises the sextic by one of its coefficients, which throws away a
// scalar; so it determines the curve only up to a quadratic twist, and the twist is what this
// script supplies.  By Eichler-Shimura the right curve has #C(F_p) = p + 1 - Tr a_p(f), so at each
// good prime the count of the untwisted model is either that or p + 1 + Tr a_p, and the pattern of
// which is the quadratic character of the twisting discriminant.  That discriminant can only ramify
// where the model is bad, so it is a signed product of those primes, and the characters give one
// linear equation over F_2 per prime.
//
// THE NEWFORM IS IDENTIFIED BY |Tr a_p|, NOT BY ATKIN-LEHNER SIGNS.  The sign convention relating
// the quaternionic and classical eigenvalues is exactly what has caused confusion here before, so
// it is avoided: the orbit is the one whose |Tr a_p| matches |p + 1 - #C(F_p)| at every good prime,
// which is convention-free and also checks that the curve is in the right isogeny class at all.
//
// ⚠ Trace(HeckeOperator(F, p)) is TWICE Tr_{K/Q}(a_p), since T_p acts on the four-dimensional
// homology as a_p on K^2.
// ⚠ Primes with Tr a_p = 0 carry no information about the character and must be dropped.

AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
D := StringToInteger(DD); N := StringToInteger(NN);
cap := StringToInteger(PRIMECAP);
PQ<u> := PolynomialRing(Rationals());
G := PQ ! (eval SEXTIC);
printf "=== level %o, sextic y^2 = %o\n", D*N, G;
error if Degree(G) notin {5, 6} or Discriminant(G) eq 0, "not a nonsingular quintic or sextic";

// the primes where this model may be bad, hence the only ones the twist can ramify at
cs := Coefficients(G);
dd := LCM([ Denominator(c) : c in cs ]);
supp := Set(PrimeDivisors(D*N))
        join Set(PrimeDivisors(Integers() ! (dd^(2*Degree(G)) * Discriminant(G))));
for c in cs do
    if c ne 0 then supp join:= Set(PrimeDivisors(Numerator(c))); end if;
    supp join:= Set(PrimeDivisors(Denominator(c)));
end for;
supp := Sort(SetToSequence(supp));
printf "the twist can only ramify at %o\n", supp;

function CountAt(poly, p)
    den := LCM([ Denominator(c) : c in Coefficients(poly) ] cat [1]);
    if den mod p eq 0 then return -1; end if;
    g := ChangeRing(PQ ! (den^2 * poly), GF(p));
    if Degree(g) lt Degree(PQ ! (den^2*poly)) - 1 or Degree(g) lt 3 or Discriminant(g) eq 0 then
        return -1;
    end if;
    return #Points(HyperellipticCurve(g));
end function;

// the counts of the untwisted model
ps := [ p : p in PrimesInInterval(3, cap) | (D*N) mod p ne 0 ];
cnt := AssociativeArray();
for p in ps do
    c := CountAt(G, p);
    if c ge 0 then cnt[p] := c; end if;
end for;
gp := Sort([ p : p in Keys(cnt) ]);
printf "%o good primes below %o\n", #gp, cap;

// the newform orbit, by |Tr a_p|
MS := ModularSymbols(D*N, 2);
NS := NewSubspace(CuspidalSubspace(MS));
hit := 0; traces := AssociativeArray();
for F in NewformDecomposition(NS) do
    if Dimension(F) ne 4 then continue; end if;
    tr := AssociativeArray(); good := true;
    for p in gp do
        t := (Integers() ! Trace(HeckeOperator(F, p))) div 2;
        tr[p] := t;
        if Abs(p + 1 - cnt[p]) ne Abs(t) then good := false; break; end if;
    end for;
    if good then hit := F; traces := tr; break; end if;
end for;
if hit cmpeq 0 then
    printf "=> no dimension-two newform orbit at level %o has |Tr a_p| = |p + 1 - #C(F_p)| at every\n", D*N;
    printf "   good prime, so this curve is not in the isogeny class of any of them\n";
    quit;
end if;
printf "matched a dimension-two newform orbit on |Tr a_p| at all %o primes\n", #gp;

// one equation over F_2 per prime with Tr a_p nonzero
rows := []; rhs := []; used := [];
for p in gp do
    t := traces[p];
    if t eq 0 then continue; end if;
    // 0 when the untwisted count already matches, 1 when the sign is opposite
    b := (p + 1 - cnt[p] eq t) select 0 else 1;
    Append(~rows, [ GF(2) | (LegendreSymbol(q, p) eq 1 select 0 else 1) : q in supp ]
                  cat [ GF(2) | (p mod 4 eq 1 select 0 else 1) ]);
    Append(~rhs, GF(2) ! b);
    Append(~used, p);
end for;
printf "%o primes carry information (Tr a_p nonzero)\n", #used;
error if #used lt #supp + 4, "too few informative primes; raise PRIMECAP";
M := Matrix(GF(2), #rows, #supp+1, &cat rows);
ok, v := IsConsistent(Transpose(M), Vector(GF(2), rhs));
if not ok then
    printf "=> the Legendre symbols are inconsistent: no single quadratic twist fits the counts,\n";
    printf "   which means the curve is NOT a twist of the modular one\n";
    quit;
end if;
error if Rank(M) ne #supp+1, "the twist is not determined uniquely; raise PRIMECAP";
sgn := (v[#supp+1] eq GF(2)!1) select -1 else 1;
t := sgn * (&*([Integers() | 1] cat [ supp[i] : i in [1..#supp] | v[i] eq GF(2)!1 ]));
printf "\ntwisting constant: %o\n", t;

C := HyperellipticCurve(t*G);
Cm := ReducedMinimalWeierstrassModel(IntegralModel(C));
printf "the model: %o\n", Cm;
printf "\nverification at every good prime: #C(F_p) against p + 1 - Tr a_p\n";
bad := 0;
for p in gp do
    c := CountAt(t*G, p);
    want := p + 1 - traces[p];
    if c ne want then bad +:= 1; printf "   p = %-5o got %-6o want %o   MISMATCH\n", p, c, want; end if;
end for;
printf "%o mismatches over %o primes\n", bad, #gp;
printf "%o\n", bad eq 0
    select "=> the twisted model matches the newform at every good prime tested"
    else   "=> the twist does not reproduce the counts";
quit;
