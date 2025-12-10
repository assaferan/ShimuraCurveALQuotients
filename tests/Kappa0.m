procedure test_Kappa0()

    printf "Testing Kappa0...";
    Q := AssociativeArray();
    O := AssociativeArray();
    L := AssociativeArray();
    SCL := AssociativeArray();
    SCL[6] := ShimuraCurveLattice(6,1);
    SCL[10] := ShimuraCurveLattice(10,1);
    
    kappa0_data := AssociativeArray();
    // verifying [Yang, Example 21, p. 24-25] and [Err, p. 850]
    kappa0_data[6] := [ <3,-4,{<2,-8>,<3,-4>}>,
                        <1,-3,{<2,-4>}>,
                        <1,-24,{ <2, -6> }>,
                        <3,-24,{ <2, -8>, <3, -4> }>,
                        <1,-163,{ <2, -4>, <3, -11>, <7, -4>, <19,-4>, <23,-4> }>,
                        <3,-163,{ <2, -40/3>, <3, -4>, <5, -4>, <11,-4>, <17,-4> }>,
                        ];

    kappa0_data[10] := [ <3,-68,{ <2, -8>,  <5, -14/3>}>,
                        <2,-68,{ <2, -6>,  <5, -6>}>
                      ];

    for D in Keys(kappa0_data) do
        for datum in kappa0_data[D] do
            m,d,log_coeffs := Explode(datum);
            assert Kappa0(m,d,SCL[D]`Q,ElementOfNorm(SCL[D]`Q,-d,SCL[D]`O,SCL[D]`basis_L)) eq LogSum(log_coeffs);
        end for;
    end for;

    printf "Done!\n";
    return;
end procedure;

test_Kappa0();