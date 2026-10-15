// W_k is the only transcendental ingredient of Theorem 3.1 and it is evaluated once per lattice
// point, so a direct Bessel call per term dominates everything else.  W_k(x) e^x is smooth and
// slowly varying, so a panelled Chebyshev approximation of it reproduces W_k to working precision
// at the cost of one polynomial evaluation.  The panels are needed because the useful range of x
// spans a factor of fifty and a single expansion over that range would need hundreds of nodes.

function ChebFit(F, a, b, n, RR)
    pi := Pi(RR);
    xs := [ (a+b)/2 + (b-a)/2 * Cos((2*j-1)*pi/(2*n)) : j in [1..n] ];
    fs := [ F(x) : x in xs ];
    ts := [ Cos((2*j-1)*pi/(2*n)) : j in [1..n] ];
    cs := [ RR | ];
    for m in [0..n-1] do
        s := RR ! 0;
        for j in [1..n] do
            s +:= fs[j] * Cos(m*(2*j-1)*pi/(2*n));
        end for;
        Append(~cs, 2*s/n);
    end for;
    return cs;
end function;

function ChebEval(cs, a, b, x)
    t := (2*x - a - b)/(b - a);
    n := #cs;
    d := Parent(x) ! 0; dd := Parent(x) ! 0;
    for m in [n..2 by -1] do
        sv := d; d := 2*t*d - dd + cs[m]; dd := sv;
    end for;
    return t*d - dd + cs[1]/2;
end function;

// Panelled approximation of x -> W_k(x) on [a,b].  The panels are GEOMETRIC, not equal width:
// W_k has a singularity at x = 0, so on an equal-width first panel [0.5, 6.5] the Chebyshev
// series converges only at the rate the Bernstein ellipse through the origin allows, which cost
// twenty digits.  Panels in a fixed ratio keep the origin equally far from every one of them.
function WkPanels(Wk, a, b, RR : ratio := 1.5, deg := 40)
    r := RR ! ratio;
    nps := Maximum(1, Ceiling(Log(b/a)/Log(r)));
    edges := [ a * r^i : i in [0..nps] ];
    coeffs := [];
    F := func< x | Wk(x) * Exp(x) >;
    for i in [1..nps] do
        Append(~coeffs, ChebFit(F, edges[i], edges[i+1], deg, RR));
    end for;
    return <edges, coeffs, Log(r), a, nps>;
end function;

function WkFast(P, x)
    edges := P[1]; coeffs := P[2]; lr := P[3]; a := P[4]; nps := P[5];
    i := Minimum(nps, Maximum(1, Ceiling(Log(x/a)/lr)));
    return ChebEval(coeffs[i], edges[i], edges[i+1], x) * Exp(-x);
end function;
