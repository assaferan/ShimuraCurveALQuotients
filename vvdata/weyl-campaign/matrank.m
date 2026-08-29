// Is the rank deficiency already in the DIVISOR MATRIX, before any weakly-holomorphic
// restriction?  coeffs_to_divisor_matrix builds a 0/1 matrix whose column for discriminant d
// is the indicator of c(d)*{squares}, where c(d) = d/4 if 4|d else d.  Two discriminants with
// the SAME CORE c(d) therefore give IDENTICAL columns -- a rank deficiency needing no eta
// quotients, no polymake, no CM points.
//
//   magma -b matrank.m
//
// Measured pole orders come from the spanprobe runs:
//   38_5 (FAILS, deficit 1): poleord 190, cols 36
//   38_7 (works, deficit 0): poleord 134 cols 25 ; poleord 266 cols 49
//   34_3 (works, deficit 0): poleord  51 cols 12 ; poleord 102 cols 22
AttachSpec("ShimuraQuotients.spec");

cases := [ <38, 5, 190>, <38, 7, 134>, <38, 7, 266>, <34, 3, 51>, <34, 3, 102> ];

for c in cases do
    D := c[1]; N := c[2]; P := c[3];
    min_m := -P;
    // num_coeffs must reach row 1 - min_m - m for m = 1 (the shallowest pole); take a margin.
    for nc in [P + 1, 2*P] do
        mat, rds := ProbeDivisorMatrix(min_m, D, N, nc);
        r := Rank(mat);
        printf "MATRANK %o_%o poleord %o num_coeffs %o nds %o cols %o rank %o deficit %o\n",
               D, N, P, nc, #rds, Ncols(mat), r, Ncols(mat) - r;
    end for;
    // Which discriminants collide under the core map c(d)?
    mat, rds := ProbeDivisorMatrix(min_m, D, N, 2*P);
    cores := AssociativeArray();
    for d in rds do
        cd := (d mod 4 eq 0) select d div 4 else d;
        if not IsDefined(cores, cd) then cores[cd] := []; end if;
        Append(~cores[cd], d);
    end for;
    coll := [<k, cores[k]> : k in Sort([x : x in Keys(cores)]) | #cores[k] gt 1];
    printf "MATCORE %o_%o poleord %o ncollisions %o %o\n", D, N, P, #coll, coll;
end for;
printf "DONE\n";
quit;
