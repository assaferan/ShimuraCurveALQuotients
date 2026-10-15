// Stage 2.  What does one evaluation by Theorem 3.1 cost, and how does that grow with the base?
//
// Nelson's formula needs NO fundamental domain and no linear algebra, so the probe runs at bases
// where Voight-Willis cannot even compute the domain.  Three things are measured per base:
//
//   * the number of lattice points summed over, which is what the running time is proportional to;
//   * the largest n at which the Fourier coefficient a_f(n) of the Jacquet-Langlands correspondent
//     is required, which is a separate and independent cost;
//   * the relative change in the value between two truncations, which says what precision the
//     cutoff actually buys.
//
//   magma -b BASES:=6_1,10_1,15_1 K:=2 PR:=30 CUT:=40 CUT0:=32 cost.m < /dev/null
SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/nelson/nelson2.m";

k := StringToInteger(K); PREC := StringToInteger(PR);
cutv := StringToInteger(CUT); cut0 := StringToInteger(CUT0);
CC := ComplexField(PREC); RR := RealField(PREC); Ic := CC.1;
bases := [ [StringToInteger(y) : y in Split(x, "_")] : x in Split(BASES, ",") ];

// a_f(n) for a weight-2 newform of level DN, taken from an elliptic curve of that conductor
// when one exists -- a_p is then a point count and 10^5 coefficients are free.  Otherwise the
// run reports the term count and timing only, with the coefficients replaced by a surrogate of
// the right size; those rows say nothing about precision and are marked.
function Coeffs(D, N, k, nmax, CC)
    DN := D*N;
    if k eq 2 then
        try
            Es := EllipticCurves(CremonaDatabase(), DN);
        catch e
            Es := [];
        end try;
        if #Es gt 0 then
            E := Es[1];
            L := LSeries(E);
            cs := LGetCoefficients(L, nmax);
            return [ CC ! cs[n] : n in [1..nmax] ], "elliptic curve of conductor " cat IntegerToString(DN);
        end if;
    end if;
    // surrogate: multiplicative, Ramanujan-size, with the forced values at p | DN
    cs := [ CC | 0 : n in [1..nmax] ];
    ap := AssociativeArray();
    for p in PrimesUpTo(nmax) do
        if DN mod p eq 0 then
            ap[p] := CC ! ((-1)^(p mod 3) * p^(k div 2 - 1));
        else
            th := 2*Pi(RealField(PREC))*Sqrt(RealField(PREC) ! p);
            ap[p] := CC ! (2 * Cos(th) * p^((k-1)/2));
        end if;
    end for;
    cs[1] := CC ! 1;
    for n in [2..nmax] do
        fa := Factorisation(n);
        v := CC ! 1;
        for pe in fa do
            p := pe[1]; e := pe[2];
            // a_{p^e} from the Hecke recursion
            prev2 := CC ! 1; prev1 := ap[p];
            if e eq 1 then w := prev1;
            else
                w := prev1;
                for j in [2..e] do
                    nw := ap[p]*prev1 - (DN mod p eq 0 select CC!0 else (CC ! p^(k-1))*prev2);
                    prev2 := prev1; prev1 := nw; w := nw;
                end for;
            end if;
            v *:= w;
        end for;
        cs[n] := v;
    end for;
    return cs, "SURROGATE coefficients (timing and counts only)";
end function;

printf "%-9o %-6o %-9o %-5o %-8o %-11o %-9o %-8o %-11o %-11o %o\n",
   "base", "DN", "s2/DN", "#div", "setup s", "lattice pts", "eval s", "n_max",
   "rel.chg", "cancel", "a_f source";
for bn in bases do
    D := bn[1]; N := bn[2]; DN := D*N;
    s2 := &+[ d^2 : d in Divisors(DN) ];
    t0 := Cputime();
    dat := NelsonSetup(D, N, k, PREC, cutv + 2);
    tsetup := Cputime(t0);
    // the a_f range follows from the cutoff: nrd(alpha) <= (cut/4pi)^2 and n = d nrd(alpha)
    C := RR ! (cutv/(4*Pi(RR)));
    nmax := Max(50, Ceiling(DN * C^2) + 10);
    cs, src := Coeffs(D, N, k, nmax, CC);
    af := func< n | n le nmax select cs[n] else 0 >;
    z1 := CC ! (0.13 + 1.07*Ic);
    z2 := CC ! (-0.21 + 0.88*Ic);
    t1 := Cputime();
    S, nt, nmx, big := NelsonEval(dat, z1, z2, af, RR!cutv, CC);
    teval := Cputime(t1);
    S0 := NelsonEval(dat, z1, z2, af, RR!cut0, CC);
    rel := Abs(S - S0)/Abs(S);
    printf "%-9o %-6o %-9o %-5o %-8o %-11o %-9o %-8o %-11o %-11o %o\n",
        IntegerToString(D) cat "_" cat IntegerToString(N), DN, RealField(5)!(s2/DN),
        #Divisors(DN), RealField(4)!tsetup, nt, RealField(4)!teval, nmx,
        RealField(3)!rel, RealField(4)!(big/Abs(S)), src;
end for;
printf "COST DONE\n";
quit;
