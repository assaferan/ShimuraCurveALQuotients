// What L_- looks like at 2 when N divides the CONDUCTOR of d (d = -12, -60 on X_0^15(2)), and the
// m = 0 local factor there.  Counting normalisation, as in lem:m0loc:
//   W_{0,2}(s, char(mu+L)) = (1-X) G(X),  alpha_k = 2^(-k) #{x in (mu+L)/2^k L : Q(x) = 0 mod 2^k}.
P<X> := PolynomialRing(Rationals());
KMAX := 7;
// the two shapes: 2*H (even, split) as at fundamental d, and 2*diag(u1,u2) with u1 + u2 = 0 mod 4
// (odd type) as measured at d = -12 and -60
shapes := [ <"2*H (fundamental d)", Matrix(Rationals(), 2, 2, [0,2, 2,0])>,
            <"2*diag(1,3)  u1+u2=0 (4)", Matrix(Rationals(), 2, 2, [2,0, 0,6])>,
            <"2*diag(1,-1) u1+u2=0 (4)", Matrix(Rationals(), 2, 2, [2,0, 0,-2])>,
            <"2*diag(1,1)  u1+u2=2 (4)", Matrix(Rationals(), 2, 2, [2,0, 0,2])> ];
for sh in shapes do
    nm := sh[1]; G := sh[2];
    gi := G^(-1);
    printf "\n=== %o,  det %o\n", nm, Integers()!Determinant(G);
    // cosets v in L^v/L: v = (a/2, b/2) works for all these (L^v = (1/2)L here)
    for ab in [<0,0>, <1,0>, <0,1>, <1,1>] do
        v := Vector(Rationals(), [ab[1]/2, ab[2]/2]);
        r := (v*G, v)/2;
        isiso := r eq Floor(r);
        alphas := [];
        for k in [1..KMAX] do
            cnt := 0;
            for u, w in [0..2^k-1] do
                x := Vector(Rationals(), [ab[1]/2 + u, ab[2]/2 + w]);
                q := (x*G, x)/2;
                if IsIntegral(q) and (Integers()!q) mod 2^k eq 0 then cnt +:= 1; end if;
            end for;
            Append(~alphas, cnt / 2^k);
        end for;
        W := (1 - X) * (1 + &+[P | alphas[k]*X^k : k in [1..KMAX]]);
        printf "   mu = (%o/2, %o/2): Q = %-6o %o   alphas %o\n      (1-X)G = %o\n",
               ab[1], ab[2], r, isiso select "ISOTROPIC" else "anisotropic", alphas[1..5], W;
    end for;
end for;
