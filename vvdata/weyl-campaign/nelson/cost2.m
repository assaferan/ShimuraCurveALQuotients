// Stage 2, large bases.  Same measurement as cost.m with two changes that it needs to reach the
// bases of the open quotients: the W_k table is built ONCE (it depends on the weight and the
// cutoff, not on the base), and the lattice points come from a ShortVectorsProcess, because one
// truncation at D N = 4830 holds several million of them and the list alone exhausts memory.
//
//   magma -b BASES:=330_1,462_1,510_1 K:=2 PR:=30 CUT:=32 CUT0:=26 NF:=3 cost2.m < /dev/null
//
// NF is how many eigenforms are evaluated in the SAME enumeration: a canonical model of a genus-g
// quotient needs g of them at each point, and they share everything but the coefficients, so this
// is the figure that decides what a model costs, not the single-form time.
SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/nelson/nelson3.m";

k := StringToInteger(K); PREC := StringToInteger(PR);
cutv := StringToInteger(CUT); cut0 := StringToInteger(CUT0);
nf := StringToInteger(NF);
CC := ComplexField(PREC); RR := RealField(PREC); Ic := CC.1;
bases := [ [StringToInteger(y) : y in Split(x, "_")] : x in Split(BASES, ",") ];

// Coefficients.  At weight 2 an elliptic curve of conductor D N gives a genuine Jacquet-Langlands
// correspondent and its a_p are point counts, so tens of thousands of them are free; that is the
// only case in which the convergence column means anything.  Otherwise a multiplicative surrogate
// of Ramanujan size stands in, which fixes the term count and the timing but not the value.
function Coeffs(D, N, k, nmax, howmany, CC)
    DN := D*N; out := []; src := "";
    if k eq 2 then
        try
            Es := EllipticCurves(CremonaDatabase(), DN);
        catch e
            Es := [];
        end try;
        for E in Es do
            if #out ge howmany then break; end if;
            cs := LGetCoefficients(LSeries(E), nmax);
            Append(~out, [ CC ! cs[n] : n in [1..nmax] ]);
        end for;
        if #out gt 0 then src := IntegerToString(#out) cat " elliptic curve(s) of conductor "
                                 cat IntegerToString(DN); end if;
    end if;
    j := 0;
    while #out lt howmany do
        j +:= 1;
        ap := AssociativeArray();
        for p in PrimesUpTo(nmax) do
            if DN mod p eq 0 then ap[p] := CC ! ((-1)^((p+j) mod 3) * p^(k div 2 - 1));
            else ap[p] := CC ! (2*Cos(2*Pi(RR)*Sqrt(RR ! (p+j))) * p^((k-1)/2)); end if;
        end for;
        cs := [ CC | 0 : n in [1..nmax] ]; cs[1] := CC ! 1;
        for n in [2..nmax] do
            v := CC ! 1;
            for pe in Factorisation(n) do
                p := pe[1]; e := pe[2];
                pr2 := CC ! 1; pr1 := ap[p]; w := pr1;
                for t in [2..e] do
                    w := ap[p]*pr1 - (DN mod p eq 0 select CC!0 else (CC ! p^(k-1))*pr2);
                    pr2 := pr1; pr1 := w;
                end for;
                v *:= w;
            end for;
            cs[n] := v;
        end for;
        Append(~out, cs);
        src cat:= (src eq "" select "" else " + ") cat "surrogate";
    end while;
    return out, src;
end function;

tW := Cputime();
Wtab := WkTable(k, cutv + 2, RR);
printf "W_%o table for cut <= %o built in %o s (base independent)\n", k, cutv+2, RR!Cputime(tW);
printf "%-9o %-6o %-9o %-8o %-12o %-10o %-8o %-11o %o\n",
   "base", "DN", "s2/DN", "setup s", "lattice pts", "eval s", "n_max", "rel.chg", "a_f source";
for bn in bases do
    D := bn[1]; N := bn[2]; DN := D*N;
    s2 := &+[ d^2 : d in Divisors(DN) ];
    t0 := Cputime();
    dat := NelsonSetupT(D, N, k, PREC, cutv + 2, Wtab);
    tsetup := Cputime(t0);
    C := RR ! (cutv/(4*Pi(RR)));
    nmax := Max(50, Ceiling(DN*C^2) + 10);
    cs, src := Coeffs(D, N, k, nmax, nf, CC);
    z1 := CC ! (0.13 + 1.07*Ic);
    z2 := CC ! (-0.21 + 0.88*Ic);
    t1 := Cputime();
    S, nt, nmx := NelsonEvalMulti(dat, z1, z2, cs, RR!cutv, CC);
    teval := Cputime(t1);
    S0 := NelsonEvalMulti(dat, z1, z2, cs, RR!cut0, CC);
    rel := Abs(S[1] - S0[1])/Abs(S[1]);
    printf "%-9o %-6o %-9o %-8o %-12o %-10o %-8o %-11o %o\n",
        IntegerToString(D) cat "_" cat IntegerToString(N), DN, RealField(5)!(s2/DN),
        RealField(4)!tsetup, nt, RealField(4)!teval, nmx, RealField(3)!rel, src;
end for;
printf "COST2 DONE\n";
quit;
