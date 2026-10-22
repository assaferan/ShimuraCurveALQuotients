// The exact decomposition of every open quotient, oldforms included.  See oldform.m for the
// mathematics; this is the same computation over the whole list, cheapest level first, with each
// curve in its own try so that one heavy level cannot cost the rest.
//
//   magma -b LIST:=nold_open61.txt oldform_all.m < /dev/null
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);

function Orbits(M)
    MS := ModularSymbols(M, 2, +1);
    NS := NewSubspace(CuspidalSubspace(MS));
    out := [];
    for A in NewformDecomposition(NS) do
        k := Dimension(A);
        sg := AssociativeArray(); okall := true;
        for p in PrimeDivisors(M) do
            T := AtkinLehnerOperator(A, p); n := Nrows(T);
            if T eq IdentityMatrix(BaseRing(T), n) then sg[p] := 1;
            elif T eq -IdentityMatrix(BaseRing(T), n) then sg[p] := -1;
            else okall := false; end if;
        end for;
        if okall then Append(~out, <k, sg>); end if;
    end for;
    return out;
end function;

procedure One(D, N, W)
    g := GenusShimuraCurveQuotient(D, N, W);
    PN := PrimeDivisors(N);
    total := 0; dims := []; lev := [];
    for Np in Divisors(N) do
        M := D * Np;
        S := [ p : p in PN | M mod p ne 0 ];
        orbs := Orbits(M);
        for o in orbs do
            k := o[1]; eps := o[2];
            rows := []; rhs := [];
            for m in W do
                if m eq 1 then continue; end if;
                known := 1;
                for p in PrimeDivisors(m) do
                    if D mod p eq 0 then known *:= -eps[p];
                    elif M mod p eq 0 then known *:= eps[p]; end if;
                end for;
                Append(~rows, [ GF(2) | (m mod p eq 0) select 1 else 0 : p in S ]);
                Append(~rhs, known eq 1 select GF(2)!0 else GF(2)!1);
            end for;
            nsol := 0;
            if IsEmpty(S) then
                nsol := forall{ i : i in [1..#rhs] | rhs[i] eq 0 } select 1 else 0;
            else
                B := Transpose(Matrix(GF(2), #rows, #S, &cat rows));
                b := Vector(GF(2), rhs);
                if IsConsistent(B, b) then nsol := 2^Dimension(Kernel(B)); end if;
            end if;
            if nsol gt 0 then
                total +:= k*nsol;
                dims cat:= [ k : i in [1..nsol] ];
                lev cat:= [ M : i in [1..nsol] ];
            end if;
        end for;
    end for;
    printf "%-9o %-3o %-7o %-16o %-26o %o\n", Sprintf("%o_%o", D, N), g, D*N,
           Sprint(Sort(dims)), Sprint(lev),
           total eq g select "AGREES" else Sprintf("DISAGREES (total %o)", total);
end procedure;

lines := Split(Read(LIST), "\n");
rows := [];
for ln in lines do
    if #ln eq 0 or ln[1] eq "#" or Position(ln, "OPEN") eq 0 then continue; end if;
    fs := Split(ln, " ");
    ps := Split(fs[1], "_");
    D := StringToInteger(ps[1]); N := StringToInteger(ps[2]);
    i1 := Position(ln, "["); i2 := Position(ln, "]");
    W := { StringToInteger(x) : x in Split(Substring(ln, i1+1, i2-i1-1), ",") };
    Append(~rows, <D*N, D, N, W>);
end for;
Sort(~rows, func< a, b | a[1] - b[1] >);
printf "%-9o %-3o %-7o %-16o %-26o %o\n", "base", "g", "D*N", "decomposition", "levels of the factors", "check";
for r in rows do
    ok := true;
    try
        One(r[2], r[3], r[4]);
    catch e
        ok := false;
        printf "%-9o  refused: %o\n", Sprintf("%o_%o", r[2], r[3]), e`Object;
    end try;
end for;
quit;
