// Galois structure of the fixed points of a residual Atkin-Lehner involution w on an open
// genus-3 quotient X with X/w elliptic.  Then w has 4 fixed points, none a Weierstrass point, and
// a hyperelliptic involution iota commutes with w and with Galois, so iota is a free
// Galois-equivariant involution of the 4-point Galois set Fix(w).  On k copies of a transitive
// Galois set G/H such an involution exists iff k is even or Aut(Q(P)/Q) = N(H)/H has even order.
// So: the number of rational fixed points is even, and for each field K with Aut(K/Q) of odd order
// the number of fixed points defined over a conjugate of K is an even multiple of [K:Q].
AttachSpec("ShimuraQuotients.spec"); SetColumns(0);
cs := eval Read(DATA);
open := [X : X in cs | not assigned X`IsSubhyp and X`g eq 3];
nr := 0;
for X in open do
    D := X`D; N := X`N; W := X`W; DN := D*N;
    full := {Integers()| d : d in Divisors(DN) | GCD(d, DN div d) eq 1};
    seen := {}; line := Sprintf("X_0(%o,%o)/%o:", D, N, Sort(Setseq(W))); verdict := "";
    for w in full diff W do
        coset := {AtkinLehnerMul(w, u, DN) : u in W};
        if coset in seen then continue; end if; Include(~seen, coset);
        Wbig := W join coset;
        gw := GenusShimuraCurveQuotient(D, N, Wbig);
        if gw ne 1 then line cat:= Sprintf(" w=%o: g(X/w)=%o;", w, gw); continue; end if;
        tot := AssociativeArray();
        for m in coset do
            e := NumFixedPointsByCMOrder(D, N, m);
            for d in Keys(e) do tot[d] := (IsDefined(tot, d) select tot[d] else 0) + e[d]; end for;
        end for;
        info := []; ok := true; nrat := 0; byfield := [* *];
        for d in Keys(tot) do
            if tot[d] eq 0 then continue; end if;
            c := tot[d] / #W;
            if not IsIntegral(c) then ok := false; info cat:= [Sprintf("d=%o: non-integral count %o", d, c)]; continue; end if;
            c := Integers()!c;
            deg := DegreeOfFieldOfDefinitionOfCMPoint(X, d);
            if deg eq 0 then ok := false; info cat:= [Sprintf("d=%o: %o fixed points but no CM point on X?", d, c)]; continue; end if;
            K := FieldOfDefinitionOfCMPoint(X, d);
            aut := deg eq 1 select 1 else #Automorphisms(K);
            info cat:= [Sprintf("d=%o: %o points, field degree %o, |Aut| %o", d, c, deg, aut)];
            if deg eq 1 then nrat +:= c; end if;
            if deg gt 1 and IsOdd(aut) then
                placed := false;
                for i in [1..#byfield] do
                    if IsIsomorphic(byfield[i][1], K) then byfield[i][2] +:= c; placed := true; break; end if;
                end for;
                if not placed then Append(~byfield, [* K, c *]); end if;
            end if;
        end for;
        line cat:= Sprintf(" w=%o (4 fixed points): %o;", w, Join(info, "; "));
        if not ok then line cat:= " [ambiguous]"; continue; end if;
        bad := IsOdd(nrat) or exists{b : b in byfield | IsOdd(b[2] div Degree(b[1]))};
        if bad then verdict := Sprintf(" => NOT hyperelliptic (w=%o: %o rational fixed points%o)", w, nrat, #byfield gt 0 select ", odd-automorphism field with " cat Sprint([b[2] : b in byfield]) cat " points" else ""); end if;
    end for;
    if verdict ne "" then nr +:= 1; end if;
    print line cat verdict;
end for;
printf "FIXGAL DONE: %o of %o ruled out\n", nr, #open;
quit;
