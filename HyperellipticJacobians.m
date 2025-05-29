path := GetCurrentDirectory();
ChangeDirectory("../hyperelliptic/");
load "Hyperelliptic3.magma";
ChangeDirectory(path);

function HypWeilPolynomialsGenus3Howe(q)
    L := eightpoints(GF(q));
    M := curves(L);
    fs := [m[1] : m in M | GCD(m[1], Derivative(m[1])) eq 1];
    wps := {Numerator(ZetaFunction(HyperellipticCurve(f))) : f in fs};
    return wps;
end function;

function IsEqualInWPS(u,v,d)
    m := #d;
    S_u := [i : i in [1..m] | u[i] ne 0];
    S_v := [i : i in [1..m] | v[i] ne 0];
    if S_u ne S_v then
        return false;
    end if;
    g, c := XGCD([d[i] : i in S_u]);
    a := &*[(v[S_u[i]]/u[S_u[i]])^c[i] : i in [1..#S_u]];
    return &and[v[i]/u[i] eq a^(d[i] div g) : i in S_u];
end function;

function wps_elts(q, d)
    F := GF(q);
    m := #d;
    ret := {};
    supports := {S : S in Subsets({1..m}) | not IsEmpty(S)};
    for supp in supports do
        S := [x : x in supp];
        _, c := XGCD([d[i] : i in S]);
        // TODO : can enumerate on one coordinate less, and compute the last one.
        FS := CartesianProduct([(i in S) select {x : x in F | x ne 0} else {F!0} : i in [1..m]]);
        for u in FS do
            if (&*[u[S[i]]^c[i] : i in [1..#S]] eq 1) then
                Include(~ret, u);
            end if;
        end for;
    end for;
    return ret;
end function;

// Slower than Howe, but should at least give the correct result
function HypWeilPolynomialsGenus3Magma(q)
    F := GF(q);
    // What we would actually want is the following
    // V := ProjectiveSpace(F, [2..7]);
    // pts := Points(V)
    // But Points does not work for weighted projective space
    // The line belwo should work, but for some reason does not yield all representatives?
    if q eq 3 then
        wts := [2,4,5,7,9,12];
    elif q gt 7 then
        wts := [2..7];
    else
        error "Not Implemented for q = %o\n", q;
    end if;
    pts := wps_elts(q,wts);
    // This is the naive approach
    // F6 := CartesianPower(F,6);
    // pts := [x : x in F6];
    Shiodas := &cat[ShiodaAlgebraicInvariants([p[i] : i in [1..6]] : ratsolve := true) : p in pts];
    good_Shiodas := {S : S in Shiodas | DiscriminantFromShiodaInvariants(S) ne 0};
    curves := &cat[TwistedHyperellipticPolynomialsFromShiodaInvariants(S) : S in good_Shiodas];
    wps := {Numerator(ZetaFunction(HyperellipticCurve(c))) : c in curves};
    return wps;
end function;