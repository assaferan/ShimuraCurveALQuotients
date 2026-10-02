// Positive control for the Galois-cover route: X_0^15(4)/<w_3,w_5> -> X_0^15(1)^* is the S_3
// Belyi cover lambda -> j, with j = 256 (1 - l + l^2)^3 / (l^2 (1-l)^2), and t_1 = M(j) for the
// Mobius map M sending j = 0, 1728, oo to the star values at discriminants -3, -12, -60.  From
// models_15_1.m: t_1(-3), t_1(-12) in {0, oo} and t_1(-15), t_1(-60) in {6, 2/27}.  Over the
// (unramified) -15 point the six lambda's must be Mobius-equivalent to Tu's six values, with the
// (3,3)-fibre {w, w-bar} over j = 0 going to his +-1/sqrt(-3).
CC := ComplexField(60); i := CC.1;
s3 := Sqrt(CC!-3); s15 := Sqrt(CC!-15);
tu := [s15/5, -s15/5, (1+s15)/8, (1-s15)/8, (-1+s15)/8, (-1-s15)/8];
tufix := [1/s3, -1/s3];
w := (1 + s3)/2;   // fixed point of the 3-cycle 1/(1-l): 1 - l + l^2 = 0, j = 0
P<l> := PolynomialRing(CC);
jpoly := func<jv | 256*(1 - l + l^2)^3 - jv * l^2 * (1 - l)^2>;
// Mobius through three points: M(0) = A, M(1728) = B, M(oo) = C (C may be oo: then M(j) = A + (B-A) j/1728)
function mobius(A, B, C)
    // M(j) = (a j + b)/(c j + d) with M(oo) = a/c = C, M(0) = b/d = A, M(1728) = B
    if C cmpeq Infinity() then return func<jv | A + (B - A)*jv/1728>, func<t | 1728*(t - A)/(B - A)>; end if;
    if A cmpeq Infinity() then // M(0) = oo: d = 0: M(j) = (a j + b)/(c j); a/c = C, (1728 a + b)/(1728 c) = B
        return func<jv | C + 1728*(B - C)/jv>, func<t | 1728*(B - C)/(t - C)>;
    end if;
    if B cmpeq Infinity() then // M(1728) = oo: M(j) = (C j - 1728 A)/(j - 1728)
        return func<jv | (C*jv - 1728*A)/(jv - 1728)>, func<t | 1728*(t - A)/(t - C)>;
    end if;
    // generic: solve for a, b, c, d with d = 1
    // M(0) = b = A; M(oo) = a/c = C; M(1728) = (1728 a + A)/(1728 c + 1) = B -> 1728 a + A = 1728 B c + B, a = C c
    // -> 1728 C c + A = 1728 B c + B -> c = (B - A)/(1728 (C - B))
    c := (B - A)/(1728*(C - B)); a := C*c; b := A;
    return func<jv | (a*jv + b)/(c*jv + 1)>, func<t | (t - b)/(a - c*t)>;
end function;
// Mobius T through three points (z_i -> w_i)
function mob3(z, wv)
    // T = S_w^-1 o S_z, with S_z the cross-ratio map z1, z2, z3 -> 0, 1, oo
    Sz := func<x | ((x - z[1])*(z[2] - z[3]))/((x - z[3])*(z[2] - z[1]))>;
    Swinv := func<y | (wv[1]*(wv[2] - wv[3]) - y*wv[3]*(wv[2] - wv[1]))/((wv[2] - wv[3]) - y*(wv[2] - wv[1]))>;
    return func<x | Abs(x - z[3]) lt 1e-40 select wv[3] else Swinv(Sz(x))>;
end function;
near := func<x, S | exists{y : y in S | Abs(x - y) lt 1e-25}>;
for assign in [<0, Infinity(), 6, 2/27>, <0, Infinity(), 2/27, 6>, <Infinity(), 0, 6, 2/27>, <Infinity(), 0, 2/27, 6>] do
    A := assign[1]; B := assign[2]; C := assign[3]; t15 := assign[4];
    M, Minv := mobius(A, B, C);
    j15 := Minv(CC!t15);
    lams := [r[1] : r in Roots(jpoly(j15))];
    // try T: w -> 1/s3, w-bar -> -1/s3, lams[1] -> each Tu value
    found := false; best := RealField(60)!1e9;
    for k in [1..6] do
        for ord in [[1,2],[2,1]] do
            T := mob3([w, ComplexConjugate(w), lams[1]], [tufix[ord[1]], tufix[ord[2]], tu[k]]);
            dev := Maximum([Minimum([Abs(T(lm) - y) : y in tu]) : lm in lams]);
            if dev lt best then best := dev; end if;
            if dev lt 1e-20 then found := true; end if;
        end for;
    end for;
    printf "  best deviation %o; lambda-six: %o\n", RealField(5)!best, [RealField(6)!Re(x) : x in lams];
    printf "t(-3)=%o t(-12)=%o t(-60)=%o t(-15)=%o : j(-15) = %o ; Tu's values reproduced: %o\n",
           A, B, C, t15, RealField(8)!Re(j15), found;
end for;
