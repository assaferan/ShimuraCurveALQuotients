// GAP 5, step 2: refine the local key, and pin the support rule.
//
// Step 1 (ctsw15.log) showed two things.  (a) a_E(m) = 0 unless (-m|p) != 1 for
// BOTH p | D -- i.e. unless every ramified prime of D is non-split in
// Q(sqrt(-m)), which is exactly the condition for that field to embed in B_D.
// (b) a_E(m)/H(4m) is constant on 20 of 27 classes keyed by
// ((-m|3),(-m|5),v2,v3,v5) but not all.
//
// The key was missing the 2-adic coordinate: 2 is the EICHLER prime at 15_2, so
// its local factor should depend on the 2-adic square class of m -- m mod 8 --
// not merely on v2.  Restrict to squarefree m (where H(4m) is a plain class
// number and the f-sums cannot confound the ratio) and key on m mod 8 as well.
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/ctlat.m";
DEPTH := 300;

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
bad := [ n : n in [1..90] | (n mod 4 in {1,2} and Coefficient(T3,n) ne 12*Hur(4*n))
                         or (n mod 8 eq 3 and Coefficient(T3,n) ne 24*Hur(n)) ];
printf "Hurwitz gate: %o failures\n", #bad;
error if #bad ne 0, "Hurwitz routine wrong";

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

v1 := thetabar(5, 2); vs := thetabar(5, 6); v3 := thetabar(30, 1);
G := [ Rationals() | -(1/2)*v1[i] + vs[i] - (1/2)*v3[i] : i in [1..DEPTH] ];

// (a) the support rule, stated and checked over the whole range
viol := [ m : m in [1..DEPTH-1] | G[m+1] ne 0 and
          (KroneckerSymbol(-m,3) eq 1 or KroneckerSymbol(-m,5) eq 1) ];
zero_allowed := [ m : m in [1..DEPTH-1] |
          KroneckerSymbol(-m,3) ne 1 and KroneckerSymbol(-m,5) ne 1 ];
printf "\nSUPPORT: a_E(m) = 0 whenever (-m|3) = 1 or (-m|5) = 1 : %o violations in m < %o\n",
    #viol, DEPTH;
printf "  of the %o admissible m, %o are nonzero\n", #zero_allowed,
    #[ m : m in zero_allowed | G[m+1] ne 0 ];

// (b) the refined local key, squarefree m only
printf "\n=== squarefree m: a_E(m)/H(4m) keyed by (m mod 8, (-m|3), (-m|5)) ===\n";
rat := AssociativeArray();
for m in [1..DEPTH-1] do
    if not IsSquarefree(m) or G[m+1] eq 0 then continue; end if;
    key := <m mod 8, KroneckerSymbol(-m,3), KroneckerSymbol(-m,5)>;
    if not IsDefined(rat, key) then rat[key] := {@ @}; end if;
    Include(~rat[key], G[m+1]/Hur(4*m));
end for;
nc := 0;
for key in Sort([ k : k in Keys(rat) ]) do
    if #rat[key] eq 1 then nc +:= 1; end if;
    printf "  m=%o mod 8, (-m|3)=%2o, (-m|5)=%2o :  %o distinct  %o\n",
        key[1], key[2], key[3], #rat[key],
        #rat[key] le 5 select [ x : x in rat[key] ] else "(many)";
end for;
printf "\n%o of %o refined classes carry a CONSTANT ratio\n", nc, #Keys(rat);
printf "DONE\n";
quit;
