// LEMMA, half 1: the LEVEL of every Gross lattice appearing in the theorem.
//
// For theta_L to be a form on Gamma_0(4DN) we need the level of the quadratic
// form Q = Nrd on the trace-zero part of an Eichler order of level R in B_Delta
// to divide 4DN, for every (Delta,R) the alternating sum uses.  The level is
// the least l > 0 with l*A^{-1} integral and EVEN on the diagonal.
//
// Also reported: det A, against the guess 2(Delta R)^2, and the local exponent
// structure, since the paper's local analysis (anisotropic plane at p | Delta,
// hyperbolic plane at p | R) should predict the odd part of the level exactly.
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/ctlat.m";

formlevel := function(A)
    B := A^(-1);
    d0 := LCM([ Integers() | Denominator(x) : x in Eltseq(B) ]);
    for mult in [1, 2] do
        l := mult*d0;
        C := l*B;
        if forall{ x : x in Eltseq(C) | Denominator(x) eq 1 }
           and forall{ i : i in [1..Nrows(C)] | IsEven(Integers()!C[i][i]) } then
            return l;
        end if;
    end for;
    error "level not found within 2*d0";
end function;

// every (Delta, R) the general theorem can call for at the base (D,N)
terms := function(D, N)
    out := {@ @};
    for s in PrimeDivisors(D) do
        Dp := D div s;
        for T in Subsets(SequenceToSet(PrimeDivisors(N))) do
            key := &*[ Integers() | p : p in T ];
            Rr := N div key;
            if IsOdd(#T) then
                Include(~out, <Dp*s*key, Rr>);
            else
                Include(~out, <Dp*key, Rr>);
                Include(~out, <Dp*key, Rr*s>);
            end if;
        end for;
    end for;
    return out;
end function;

bases := [ <15,2>, <21,2>, <22,3>, <33,2>, <35,2>, <55,2>, <77,2>, <10,7>,
           <34,3>, <6,5>, <26,3>, <39,2>, <35,6>, <15,14>,
           <15,1>, <210,1>, <330,1> ];

bad := 0; tot := 0;
for b in bases do
    D := b[1]; N := b[2]; M := 4*D*N;
    printf "\n=== D=%o N=%o   4DN = %o ===\n", D, N, M;
    for t in terms(D, N) do
        Dl := t[1]; Rl := t[2];
        if #PrimeDivisors(Dl) mod 2 eq 0 then continue; end if;   // indefinite: no Gross lattice
        A := ChangeRing(grossGram(Dl, Rl), Rationals());
        AZ := ChangeRing(A, Integers());
        lev := formlevel(A);
        dt := Determinant(AZ);
        tot +:= 1;
        okl := M mod lev eq 0;
        if not okl then bad +:= 1; end if;
        printf "  Gr(%o,%o): det = %-10o (2(DR)^2 = %-10o %o)  level = %-8o %o\n",
            Dl, Rl, dt, 2*(Dl*Rl)^2,
            dt eq 2*(Dl*Rl)^2 select "ok" else "DIFFERS",
            lev, okl select Sprintf("| 4DN  (4DN/lev = %o)", M div lev)
                       else "*** DOES NOT DIVIDE 4DN ***";
    end for;
end for;
printf "\n%o of %o lattices have level dividing 4DN\n", tot - bad, tot;
printf "DONE\n";
quit;
