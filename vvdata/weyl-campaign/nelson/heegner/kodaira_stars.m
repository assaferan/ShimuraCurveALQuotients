// Kodaira types of genus-one star curves from the metric dual graph of the Cerednik-Drinfeld model at
// the primes of D: the minimal metric graph of a genus-one quotient is a cycle, and the type is I_n with
// n the total length. Calibrated on stars whose curve is known (142_1 = 142a1, 158_1 = 158a1, 38_5 = 190a1),
// then applied to 14_17 (238a1 or 238a2) and 94_3 (282b1 or 282b2).
AttachSpec("ShimuraQuotients.spec"); SetColumns(0);
cases := [ <142, 1, {1,2,71,142}, ["142a1"]>, <158, 1, {1,2,79,158}, ["158a1"]>, <38, 5, {1,2,19,38,5,10,95,190}, ["190a1"]>,
           <14, 17, {1,2,7,14,17,34,119,238}, ["238a1", "238a2"]>, <94, 3, {1,2,47,94,3,6,141,282}, ["282b1", "282b2"]> ];
for c in cases do
    D, N, W, labs := Explode(c); W := {Integers() | w : w in W};
    printf "\n=== X_0^%o(%o)/W, genus %o ===\n", D, N, GenusShimuraCurveQuotient(D, N, W);
    for p in PrimeDivisors(D) do
        data := DualGraphData(D, N, p);
        nv, edges, lens, resid := DualGraphMetricQuotientFromData(data, p, W);
        nv2, e2, l2, rep := MetricGraphModel(nv, edges, lens);
        printf "p = %o: quotient graph %o vertices %o edges, minimal model %o vertices %o edges, lengths %o, total %o\n", p, nv, #edges, nv2, #e2, l2, &+l2;
    end for;
    for l in labs do
        E := EllipticCurve(l);
        printf "%o: %o\n", l, [<p, KodairaSymbol(E, p)> : p in PrimeDivisors(D)];
    end for;
end for;
