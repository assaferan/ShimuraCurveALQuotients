procedure test_genera()
//Verify that the genera are being computed correctly by comparing to the data file
    r := GetLargestPrimeIndex();
    assert r eq 7;
    star_curves := FindPairs(r); // time : 1.980
    assert #star_curves eq 2342;

    UpdateGenera(~star_curves); // time: 12
    UpdateByGenus(~star_curves);
    FilterByTrace(~star_curves); // time : 2850.270
    HHProposition1(~star_curves);
    curves := GetQuotientsAndGenera(star_curves);
    read_curves := eval Read("data/all_curves_progress.dat");
    assert &and[(<curves[j]`D, curves[j]`N, curves[j]`W> eq <read_curves[j]`D, read_curves[j]`N, read_curves[j]`W>) : j in [1..#curves]];
    assert &and[(curves[j]`g eq read_curves[j]`g) : j in [1..#curves]];
end procedure;

procedure test_data_file()
//Verify that the genera are being computed correctly by comparing to the data file
    r := GetLargestPrimeIndex();
    assert r eq 7;
    star_curves := FindPairs(r); // time : 1.980
    assert #star_curves eq 2342;
    UpdateGenera(~star_curves); // time: 12
    UpdateByGenus(~star_curves);
    FilterByTrace(~star_curves); // time : 2850.270
    HHProposition1(~star_curves);
    curves := GetQuotientsAndGenera(star_curves);
    UpdateByGenus(~curves);
    UpwardClosure(~curves);
    DownwardClosure(~curves);
    FilterByALFixedPointsOnQuotient(~curves);
    UpwardClosure(~curves);
    Genus3CoversGenus2(~curves);
    DownwardClosure(~curves);
    FilterByComplicatedALFixedPointsOnQuotient(~curves); // 64.840
    UpwardClosure(~curves);
    FilterByTrace(~curves);
    UpwardClosure(~curves);
    UpdateByIsomorphisms(~curves); // 13.750
    UpwardClosure(~curves);
    DownwardClosure(~curves);
    FilterByWSPoints(~curves);
    UpdateByIsomorphisms(~curves);
    UpwardClosure(~curves);
    DownwardClosure(~curves);
    unknownstar := [ c : c in curves | IsStarCurve(c) and not assigned c`IsSubhyp];
    for p in PrimesUpTo(10) do
        FilterStarCurvesByFpAutomorphisms(unknownstar, ~curves, p, 20 );
    end for;
    UpdateByIsomorphisms(~curves);
    UpwardClosure(~curves);
    DownwardClosure(~curves);
    FilterByWeilPolynomial(~curves : genera := {3,4,5});
    UpdateByIsomorphisms(~curves);
    UpwardClosure(~curves);
    DownwardClosure(~curves);
    FilterByDegeneracyMorphism(~curves);
    UpdateByIsomorphisms(~curves);
    UpwardClosure(~curves);
    DownwardClosure(~curves);
    read_curves := eval Read("data/all_curves_progress.dat");
    assert &and[(curves[j] eq read_curves[j]) : j in [1..#curves]];
end procedure;

