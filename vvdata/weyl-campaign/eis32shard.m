// eis32shard: EMAT rows [EFROM..ETO] of the eis32 E-constant table, full precision.
// Setup skips BorcherdsForms entirely (panel/keys/kernel not needed for the table).
//   magma -b DD:=210 NN:=1 EF:=<pool> EFROM:=1 ETO:=23 eis32shard.m
AttachSpec("ShimuraQuotients.spec");

D := 15; N := 2; keystr := "-2,-1,9,10,11,12,13,14,15";
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;
if assigned KEYS then keystr := KEYS; end if;
keys := [StringToInteger(s) : s in Split(keystr, ",")];
MAXE := 80; if assigned ME then MAXE := StringToInteger(ME); end if;

PREC := 80;
if assigned PR then PREC := StringToInteger(PR); end if;   // PR:=120 at M >= 660
CC := ComplexField(PREC); ii := CC.1; pi := Pi(CC);
ee := func< z | Exp(2*pi*ii*z) >;


M := IsOdd(D*N) select 4*D*N else 2*D*N;
ds := Divisors(M); nd := #ds;
printf "SHARD BASE %o %o  M = %o  PREC = %o\n", D, N, M, PREC;
// O(M^3) -> O(psi(M)*phi(M)); identical output, see fastcosets_check.m
load "vvdata/weyl-campaign/fastcosets.m";
reps := fastCosetReps(M);
words := [ VVSTWord(g) : g in reps ];
printf "#cosets = %o\n", #reps;

cs := Divisors(M);
ordvec := function(r)  // 24*ord at cusp c, up to positive width factors
    return [ &+[ Rationals() | r[i]*GCD(c, ds[i])^2/ds[i] : i in [1..nd] ] : c in cs ];
end function;

Epool := {@ @};
// optional external pool (python-enumerated holomorphic quotients): EF := path
// to a file containing a Magma sequence literal [ [r_1,...,r_nd], ... ]
if assigned EF then
    EXTRAE := eval Read(EF);
    for r in EXTRAE do
        rr := [ Integers() | x : x in r ];
        error if #rr ne nd, "EF vector length mismatch";
        error if exists{ o : o in ordvec(rr) | o lt 0 }, "EF vector not holomorphic";
        Include(~Epool, rr);
    end for;
    printf "#E's after external pool = %o\n", #Epool;
end if;

EFROMi := StringToInteger(EFROM); ETOi := StringToInteger(ETO);
// ---- cusp6 coset machinery, weight-general ------------------------------
triang := function(g, d)
    g2 := Matrix(Integers(), 2, 2, [d*g[1][1], d*g[1][2], g[2][1], g[2][2]]);
    c1 := g2[1][1]; c2 := g2[2][1];
    h := GCD(c1, c2);
    p1 := c1 div h; p2 := c2 div h;
    gg, u, v := XGCD(p1, p2);
    error if gg ne 1, "not primitive";
    gd := Matrix(Integers(), 2, 2, [p1, -v, p2, u]);
    sd := gd^(-1) * g2;
    a := sd[1][1]; b := sd[1][2]; e := sd[2][2];
    if a lt 0 then a := -a; b := -b; e := -e; end if;
    if e lt 0 then gd := -gd; sd := -sd; a := sd[1][1]; b := sd[1][2]; e := sd[2][2]; end if;
    return a, b, e, gd;
end function;

slashdata := function(word, tau)
    z := tau; factor := CC!1;
    for i := #word to 1 by -1 do
        if word[i][1] eq "S" then factor /:= Sqrt(z); z := -1/z;
        else z := z + word[i][2]; end if;
    end for;
    return factor, z;
end function;

tau0 := CC!0.31 + CC!1.31*ii;
tau1 := CC!(-0.57) + CC!1.73*ii;

SS := PowerSeriesRing(CC); t := SS.1;

// a0 of the eta monomial with exponent vector r (weight (sum r)/2) at coset w;
// if depth0 > 0 also print the first expansion coefficients (for E* readout).
a0at := function(w, r, wt)
    g := VVWordMatrix(w);
    tri := [ ];
    for d in ds do
        a, b, e, gd := triang(g, d);
        Append(~tri, <a, b, e, gd>);
    end for;
    W := LCM([ tri[i][3] : i in [1..nd] ]);
    L := &+[ Integers() | r[i]*tri[i][1]*(W div tri[i][3]) : i in [1..nd] ];
    if L gt 0 then return CC!0, L, W; end if;
    depth := -L + 1;
    units := [ ];
    for i->d in ds do
        a, b, e := Explode([tri[i][1], tri[i][2], tri[i][3]]);
        step := 24*a*(W div e);
        u := SS!1 + O(t^depth);
        n := 1;
        while n*step lt depth do
            u *:= 1 - ee(CC!(n*b/e))*t^(n*step);
            n +:= 1;
        end while;
        Append(~units, u);
    end for;
    produ := &*[ SS | units[i]^(r[i]) : i in [1..nd] | r[i] ne 0 ];
    c0 := Coefficient(produ, -L);
    if c0 eq 0 then return CC!0, L, W; end if;
    fac0, z0 := slashdata(w, tau0);
    fac1, z1 := slashdata(w, tau1);
    num0 := fac0^wt * &*[ CC | DedekindEta(d*z0)^(r[i]) : i->d in ds | r[i] ne 0 ];
    num1 := fac1^wt * &*[ CC | DedekindEta(d*z1)^(r[i]) : i->d in ds | r[i] ne 0 ];
    sfun := func< tau | ee(tau*L/(24*W)) *
        &*[ CC | ( DedekindEta((tri[i][1]*tau + tri[i][2])/tri[i][3]) *
                   ee(-(tri[i][1]*tau + tri[i][2])/(24*tri[i][3])) )^(r[i])
            : i in [1..nd] | r[i] ne 0 ] >;
    k0 := num0 / sfun(tau0);
    k1 := num1 / sfun(tau1);
    // RELATIVE consistency check, as in eis32k/eis32s: at M >= 660 the eta
    // evaluations lose enough digits that the old ABSOLUTE 1e-30 bound fires
    // spuriously (dev 3.2e-30 on |k| of order 1).
    kscale := Abs(k0); if kscale lt 1 then kscale := Parent(kscale)!1; end if;
    error if Abs(k0 - k1) gt 10^(-25)*kscale, "kappa not constant at wt", wt, Abs(k0-k1);
    return k0 * c0, L, W;
end function;


nw := #words;
// Emit EMAT in EXACTLY the eis32s format -- "EMAT <E> <coset> <class> re im" at
// 50 digits -- so shard logs concatenate with a rho-only eis32s log into
// something kernelrat/eisrat consume directly, with no combiner.
R50 := RealField(50);
for Ei in [EFROMi..Minimum(ETOi, #Epool)] do
    rE := Epool[Ei];
    for wi->w in words do
        v, L, W := a0at(w, rE, 3);
        if Abs(v) lt 10^(-40) then continue; end if;
        g := VVWordMatrix(w);
        printf "EMAT %o %o %o %o %o\n", Ei, wi, GCD(g[2][1] mod M, M), R50!Re(v), R50!Im(v);
    end for;
    printf "SHARDDONE %o\n", Ei;
end for;
printf "SHARD COMPLETE %o..%o\n", EFROMi, ETOi;
quit;
