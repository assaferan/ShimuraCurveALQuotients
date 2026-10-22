// THE EXACT DECOMPOSITION OF Jac(X_0^D(N)/W), OLDFORMS INCLUDED, from the degeneracy maps.
//
//   magma -b DD:=77 NN:=10 WS:=1,2,5,7,10,11,14,22,35,55,70,77,110,154,385,770 oldform.m
//
// By Jacquet-Langlands the differentials on X_0^D(N) are S_2(DN) new at every p | D, so they
// decompose over the newform orbits of every level M with D | M | DN -- not only M = DN.  An orbit
// of level M appears through the degeneracy maps, one copy per divisor of DN/M.
//
// WHAT THE DEGENERACY MAPS DO TO THE ATKIN-LEHNER SIGNS.  For p not dividing M, the two maps into
// level pM span a w_p-stable plane on which w_p SWAPS them; its trace is therefore 0 and its
// eigenvalues are +1 and -1, one each.  So for S = the primes of N at which the orbit is OLD, the
// 2^#S copies split into simultaneous w_p-eigenspaces indexed by a free sign eta_p for each p in S,
// each sign vector occurring exactly once.  At the primes p | M the sign is the orbit's own
// eps_p, and at p | D it is MINUS that, the Jacquet-Langlands flip.
//
// So the orbit contributes dim(f) times the number of sign vectors eta making every w_m with
// m in W act as +1:
//     prod_{p in S, p | m} eta_p  *  prod_{p | m, p | D} (-eps_p)  *  prod_{p | m, p | N, p | M} eps_p  =  +1
// which is a linear system over F_2 with 0 or 2^nullity solutions.
//
// THE CHECK, and it is a real one: the total must equal the genus of the quotient, computed
// independently from the trace formula.  Nothing here is fitted to that number.

AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
D := StringToInteger(DD); N := StringToInteger(NN);
W := { StringToInteger(x) : x in Split(WS, ",") };
require IsSquarefree(D*N) : "this implementation assumes D*N squarefree (the AL operators are at exact primes)";

X := CreateShimuraQuot(D, N, W);
g := GenusShimuraCurveQuotient(D, N, W);
printf "=== X_0^%o(%o)/%o : genus %o (trace formula)\n", D, N, Sort(Setseq(W)), g;

PD := PrimeDivisors(D); PN := PrimeDivisors(N);

// newform orbits of a level, with their dimensions and Atkin-Lehner signs
function Orbits(M)
    MS := ModularSymbols(M, 2, +1);
    NS := NewSubspace(CuspidalSubspace(MS));
    out := [];
    for A in NewformDecomposition(NS) do
        k := Dimension(A);
        sg := AssociativeArray();
        okall := true;
        for p in PrimeDivisors(M) do
            T := AtkinLehnerOperator(A, p);
            n := Nrows(T);
            if T eq IdentityMatrix(BaseRing(T), n) then sg[p] := 1;
            elif T eq -IdentityMatrix(BaseRing(T), n) then sg[p] := -1;
            else okall := false; end if;
        end for;
        if okall then Append(~out, <k, sg>); else
            printf "  ! level %o: an orbit has a non-scalar w_p, skipped\n", M;
        end if;
    end for;
    return out;
end function;

total := 0;
pieces := [];
// M = D * N', N' | N
for Np in Divisors(N) do
    M := D * Np;
    S := [ p : p in PN | M mod p ne 0 ];     // primes of N where the orbit is OLD
    orbs := Orbits(M);
    // printed because a level that silently produced nothing would look exactly like a level whose
    // orbits all fail the sign conditions, and the totals could still come out right by luck
    printf "  level %-6o old at %-10o %o newform orbits, dimensions %o\n",
           M, Sprint(S), #orbs, [ o[1] : o in orbs ];
    for o in orbs do
        k := o[1]; eps := o[2];
        // build the F_2 system over the unknowns x_p, p in S
        rows := []; rhs := [];
        for m in W do
            if m eq 1 then continue; end if;
            known := 1;
            for p in PrimeDivisors(m) do
                if D mod p eq 0 then known *:= -eps[p];
                elif M mod p eq 0 then known *:= eps[p];
                end if;
            end for;
            Append(~rows, [ GF(2) | (m mod p eq 0) select 1 else 0 : p in S ]);
            Append(~rhs, known eq 1 select GF(2)!0 else GF(2)!1);
        end for;
        nsol := 0;
        if IsEmpty(S) then
            nsol := forall{ i : i in [1..#rhs] | rhs[i] eq 0 } select 1 else 0;
        else
            // Magma solves x*B = b, so the system A*x = b is posed with B = Transpose(A); the
            // solutions then differ by Kernel(B), a subspace of F_2^S.
            B := Transpose(Matrix(GF(2), #rows, #S, &cat rows));
            b := Vector(GF(2), rhs);
            solvable := IsConsistent(B, b);
            nsol := solvable select 2^Dimension(Kernel(B)) else 0;
        end if;
        if nsol gt 0 then
            total +:= k*nsol;
            Append(~pieces, <M, k, nsol, [ <p, eps[p]> : p in Sort(PrimeDivisors(M)) ]>);
        end if;
    end for;
end for;

printf "\n%-8o %-5o %-7o %o\n", "level", "dim", "copies", "Atkin-Lehner signs of the orbit";
for t in pieces do
    printf "%-8o %-5o %-7o %o\n", t[1], t[2], t[3], t[4];
end for;
printf "\nTOTAL DIMENSION %o against genus %o: %o\n", total, g,
       total eq g select "AGREES" else "DISAGREES -- the bookkeeping is wrong";
printf "decomposition (dimensions of the simple factors, with multiplicity): %o\n",
       Sort(&cat[ [ t[2] : i in [1..t[3]] ] : t in pieces ]);
quit;
