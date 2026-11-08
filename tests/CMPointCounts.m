// The weights and point counts of CM points that the dual test and the component groups use, against
// independent sources.
//   * CMPointWeight(D, N, d), the order of the stabiliser of a point of discriminant d on the star (the
//     Heegner series weights P_d by 1/e_d): units of the CM order, doubled at a fixed point of an
//     Atkin-Lehner involution.  Whether P_d is such a fixed point is read off Ogg's fixed-point counts
//     (NumFixedPointsByCMOrder), an independent routine; in particular P_-4 is fixed by w_2 (1 + i).
//   * CMSpecialisation, which places the star points of discriminant d on the components of the special
//     fibre: the number of star points it finds must be the degree of their field of definition
//     (DegreeOfFieldOfDefinitionOfCMPoint, class field theory).  This failed at P_-3 on X_0^6(247)*
//     before the units of the CM order were taken out of the stabiliser (2026-11-08).

printf "CM point weights and counts against Ogg and class field theory...";
// weights at fixed and non-fixed points: P_-4 with w_2 present, P_-3 with w_3, and three that are not doubled or are
for c in [<141, 2, -4>, <6, 247, -3>, <6, 247, -52>, <21, 34, -4>, <15, 38, -12>] do
    D, N, d := Explode(c);
    fixed := exists{m : m in Divisors(D*N) | m gt 1 and GCD(m, (D*N) div m) eq 1 and
                    IsDefined(NumFixedPointsByCMOrder(D, N, m), d) and NumFixedPointsByCMOrder(D, N, m)[d] ne 0};
    units := d eq -3 select 3 else (d eq -4 select 2 else 1);
    assert CMPointWeight(D, N, d) eq units * (fixed select 2 else 1);
end for;
assert CMPointWeight(141, 2, -4) eq 4 and CMPointWeight(6, 247, -3) eq 6;
// star points on the special fibre against the degree of the field of definition
for c in [<6, 247, [2, 3], [-3, -52, -456]>, <141, 2, [3], [-4, -12, -1128]>] do
    D, N, ps, ds := Explode(c);
    full := {Integers() | m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1};
    X := CreateShimuraQuot(D, N, full); X`g := GenusShimuraCurveQuotient(D, N, full); X`CurveID := 0;
    star := rep{s : s in GetQuotientsAndGenera([X]) | IsStarCurve(s)};
    for p in ps do
        G := ComponentGroupData(D, N, full, p);
        for d in ds do
            _, total := CMSpecialisation(G, d);
            assert total eq DegreeOfFieldOfDefinitionOfCMPoint(star, d);
        end for;
    end for;
end for;
printf " ok\n";
