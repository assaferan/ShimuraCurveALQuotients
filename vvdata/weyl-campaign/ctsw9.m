// GAP 5c: N = 1, D odd -- the six model bases the odd law can actually reach.
//
// Prop 9.15 assumes N an odd PRIME; at N = 1 its Eichler factor
// 12(N-chi_N)/(N^2-1) is 0/0, so N = 1 needs its own statement.  The general
// theorem at N = 1 has S empty and no Eichler prime at all:
//     E = Tel_s(D'; 1) = [((s+1)/2) thetabar(D',s) - thetabar(D',1)] / ((s-1)/2)
// so the expected shape is the ramified product and the conductor corrections
// with NO Eichler factor, times a constant.
//
// These are the real thing: 15_1, 35_1, 39_1, 55_1, 65_1, 87_1 are exactly the
// six bases in data/models with DN odd.  Both choices of s are run, since the
// theorem says the answer does not depend on which prime is s.
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/ctlat.m";
DEPTH := 250;

Hur := function(n)
    if n eq 0 then return Rationals()!(-1)/12; end if;
    if n mod 4 in {1, 2} then return Rationals()!0; end if;
    s := Rationals()!0; f := 1;
    while f^2 le n do
        if n mod f^2 eq 0 then
            d := -(n div f^2);
            if d mod 4 in {0, 1} then
                w := d eq -3 select 6 else (d eq -4 select 4 else 2);
                s +:= Rationals()!(2*ClassNumber(QuadraticForms(d)))/w;
            end if;
        end if;
        f +:= 1;
    end while;
    return s;
end function;
T3 := ThetaSeries(StandardLattice(3), 4*DEPTH);
error if exists{ n : n in [1..90] |
    (n mod 4 in {1,2} and Coefficient(T3,n) ne 12*Hur(4*n)) or
    (n mod 8 eq 3 and Coefficient(T3,n) ne 24*Hur(n)) }, "Hurwitz routine wrong";
printf "Hurwitz gate: passed\n";

sig := function(p, j)
    if j lt 0 then return Rationals()!0; end if;
    return Rationals()!((p^(j+1) - 1) div (p - 1));
end function;
thetabar := function(Dp, Rl)
    A := ChangeRing(grossGram(Dp, Rl), Rationals());
    L := LatticeWithGram(A);
    reps := Representatives(Genus(L));
    raw := [ Rationals() | 0 : m in [0..2*DEPTH] ];
    den := Rationals()!0;
    for Lr in reps do
        w := Rationals()!1/#AutomorphismGroup(Lr);
        T := ThetaSeries(Lr, 2*DEPTH);
        for m in [0..2*DEPTH] do raw[m+1] +:= w*Coefficient(T, m); end for;
        den +:= w;
    end for;
    raw := [ x/den : x in raw ];
    return [ raw[2*m+1] : m in [0..DEPTH-1] ];
end function;

Ds := [ 15, 35, 39, 55, 65, 87 ];   // every DN-odd base in data/models
for D in Ds do
    printf "\n=== D = %o = %o,  N = 1 ===\n", D, PrimeDivisors(D);
    for s in PrimeDivisors(D) do
        Dp := D div s;
        t1 := thetabar(Dp, 1); ts := thetabar(Dp, s);
        // Tel_s(D';1), the mass factors reduce to (s+1)/2 against 1
        E := [ Rationals() | (((s+1)/2)*ts[i] - t1[i])/((s-1)/2) : i in [1..DEPTH] ];
        R := AssociativeArray();
        for m in [1..DEPTH-1] do
            d0 := FundamentalDiscriminant(-4*m);
            c := Isqrt((-4*m) div d0);
            L := Rationals()!1;
            for p in PrimeDivisors(D) do
                k := Valuation(c,p); chi := KroneckerSymbol(d0,p);
                L *:= (1 - chi)/((p-1)*(sig(p,k) - chi*sig(p,k-1)));
            end for;
            if L eq 0 then
                if E[m+1] ne 0 then printf "   ** m=%o: L=0 but a_E=%o\n", m, E[m+1]; end if;
                continue;
            end if;
            key := <m mod 8>;
            if not IsDefined(R, key) then R[key] := {@ @}; end if;
            Include(~R[key], -E[m+1]/(Hur(4*m)*L));
        end for;
        vals := {@ @};
        for k in Keys(R) do for x in R[k] do Include(~vals, x); end for; end for;
        printf "  s = %-3o D' = %-3o :  residual constant?  %o distinct value%o  %o\n",
            s, Dp, #vals, #vals eq 1 select "" else "s",
            #vals le 6 select [ x : x in vals ] else "(many)";
        if #vals gt 1 then
            for k in Sort([ x : x in Keys(R) ]) do
                printf "        m = %o mod 8 : %o\n", k[1],
                    #R[k] le 4 select [x : x in R[k]] else "(many)";
            end for;
        end if;
    end for;
end for;
printf "DONE\n";
quit;
