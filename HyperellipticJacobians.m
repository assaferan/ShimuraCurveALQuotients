load "../hyperelliptic/Hyperelliptic3.magma";

function HypWeilPolynomialsGenus3(q)
    L := eightpoints(GF(q));
    M := curves(L);
    fs := [m[1] : m in M | GCD(m[1], Derivative(m[1])) eq 1];
    wps := {Numerator(ZetaFunction(HyperellipticCurve(f))) : f in fs};
    return wps;
end function;