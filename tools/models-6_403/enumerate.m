// GENERIC FIBRE ENUMERATION for a genus-four cover Y = X_0^D(N)/H of a genus-two star C (candidate
// model y^2 + h y = f).  For every rational Weierstrass point e of C (each gives a rational theta
// characteristic, hence a normal form), every fibre x0 of height <= BOUND and every twist c supported
// on S u {-1}: does Y_c : y^2 = q1(t^2/c), q1 the quintic with the fibre at u = 0 and e at u = oo,
// have the point counts of Y?  Uses #Y_c(F_p) = #C(F_p) + #C'_c(F_p) - (p+1), C'_c : v^2 = c u q1(u).
//   magma -b DD:=6 NN:=403 FC:="[9,-21,33,-47,15,9,1]" HC:="[0,1,1]" HH:="{1,2,13,26,93,186,1209,2418}" BOUND:=1000000 generic_enum.m
AttachSpec("ShimuraQuotients.spec"); SetColumns(0);
D := StringToInteger(DD); N := StringToInteger(NN); BND := assigned BOUND select StringToInteger(BOUND) else 1000000;
P<x> := PolynomialRing(Rationals());
C := HyperellipticCurve(P!eval FC, P!eval HC); Cs := SimplifiedModel(C); g := HyperellipticPolynomials(Cs);
H := eval HH; S := Set(PrimeDivisors(D*N));
G := { m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1 };
Yq := CreateShimuraQuot(D, N, H); Yq`g := GenusShimuraCurveQuotient(D, N, H);
Xq := CreateShimuraQuot(D, N, G); Xq`g := GenusShimuraCurveQuotient(D, N, G);
primes := [ p : p in PrimesUpTo(47) | (D*N) mod p ne 0 ];
tY := [ ComputePointsViaTrace(Yq, p, 1) : p in primes ];
mX := [ #Points(ChangeRing(Cs, GF(p))) : p in primes ];
printf "star model vs trace at %o: %o\n", primes, mX eq [ ComputePointsViaTrace(Xq, p, 1) : p in primes ];
need := [ tY[i] - mX[i] + primes[i] + 1 : i in [1..#primes] ];
pts := Points(Cs : Bound := BND);
Wf := Sort(Setseq({ pt[1]/pt[3] : pt in pts | pt[3] ne 0 and pt[2] eq 0 }));        // finite rational Weierstrass points
es := [* e : e in Wf *]; if Degree(g) eq 5 then Append(~es, Infinity()); end if;
fib := Sort(Setseq({ pt[1]/pt[3] : pt in pts | pt[3] ne 0 and pt[2] ne 0 }));
x0s := [* x0 : x0 in fib *]; if Degree(g) eq 6 then Append(~x0s, Infinity()); end if;
printf "%o rational points (height <= %o); Weierstrass points %o; %o non-Weierstrass fibres; Jac[2](Q) order %o\n",
       #pts, BND, es, #x0s, #TwoTorsionSubgroup(Jacobian(Cs));
PT<t> := PolynomialRing(Rationals());
twists := Sort([ s * &*[ Integers() | q : q in T ] : T in Subsets(S), s in [1, -1] ]);
hits := [];
for e in es do for x0 in x0s do
    // u = (x - x0)/(x - e); q1(u) = (numerator)(u) of degree 5
    if Type(e) eq Infty then                       // u = x - x0
        q1 := Evaluate(g, t + x0);
    elif Type(x0) eq Infty then                    // u = 1/(x - e)
        q1 := &+[ Coefficient(g, i) * (e*t + 1)^i * t^(6-i) : i in [0..6] ];
    else
        q1 := &+[ Coefficient(g, i) * (e*t - x0)^i * (t - 1)^(6-i) : i in [0..6] ];
    end if;
    q1 := PT ! q1;
    assert Degree(q1) eq 5 and Coefficient(q1, 0) ne 0;
    for c in twists do
        Cp := ReducedMinimalWeierstrassModel(HyperellipticCurve(c * t * q1));
        ok := true;
        for i in [1..#primes] do
            try n := #Points(ChangeRing(Cp, GF(primes[i]))); catch err n := need[i]; end try;
            if n ne need[i] then ok := false; break; end if;
        end for;
        if ok then
            Append(~hits, <e, x0, c>);
            Y := Evaluate(q1, t^2/c);
            // the primitive integral multiple of q1(t^2/c); the scaling is by a constant, so y^2 = (that)
            // is a quadratic twist of Y at worst, and verify.m settles the isomorphism over Q
            den := LCM([Denominator(a) : a in Coefficients(Y)]); Yi := PT ! (den * Y);
            cont := GCD([Integers() | a : a in Coefficients(Yi)]);
            printf "MATCH: base e = %o, fibre x0 = %o, twist c = %o\n   C' = %o\n   Y : y^2 = q1(t^2/c), primitive integral multiple %o\n", e, x0, c, Cp, Yi div cont;
        end if;
    end for;
end for; end for;
printf "%o matching (base, fibre, twist) triple(s) over %o bases, %o fibres, %o twists\n", #hits, #es, #x0s, #twists;
quit;
