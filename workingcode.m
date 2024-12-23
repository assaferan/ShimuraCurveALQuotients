// AttachSpec("shimuraquots.spec");
intrinsic GetHyperellipticCandidates(:recompute_data:=false, read_data :=true) -> SeqEnum
{.}
    SetDebugOnError(true);
    SetVerbose("ShimuraQuotients", 3);

    // Find the largest prime we need to consider for the
    // inequality in Proposition 1.

    if recompute_data then 
        r := GetLargestPrimeIndex();

        assert r eq 7;

        // Find all pairs (D,N) satisfying the inequality of
        // Proposition 1.
        time star_curves := FindPairs(r); // time : 1.980

        // I added some code that just
        // focuses on the star quotients X_0^*(D,N)

        assert #star_curves eq 2342;

        time UpdateGenera(~star_curves); // time: 12

        VerifyHHTable1(star_curves);

        UpdateByGenus(~star_curves);


        FilterByTrace(~star_curves); // time :

        VerifyHHTable2(star_curves);

        // Create a list of all Atkin-Lehner quotients
        // compute their genera, and store the covering structure.

        // writing to a file, in case we would like to load it directly

        Write("star_curves_point_count.dat", "" : Overwrite);

        for curve in star_curves do
            Write("star_curves_point_count.dat", Sprint(curve, "Magma"));
        end for;

        // testing that reading the file works
        if false then //check data
            read_curves := eval Read("star_curves_point_count.dat");
            assert #read_curves eq #star_curves;
            assert &and[(read_curves[j] eq star_curves[j]) : j in [1..#star_curves]];
            HHProposition1(~star_curves);
            VerifyHHProposition1(star_curves);
        end if;

        time curves := GetQuotientsAndGenera(star_curves); // 148.660

        // updating classification from the genera we computed
        UpdateByGenus(~curves);

        // upward closure - if X-->Y is dominant then gon(X) >= gon(Y). [Poonen, A.1.(vii)]
        // In particular, if Y is non-hyperelliptic, so is X.
        UpwardClosure(~curves);

        // downward closure - if covered by subhyperelliptic, then subhyperelliptic
        DownwardClosure(~curves);

        // Using the fact that if w acts non-trivially ans has more than
        // 4 fixed points on X, then X is non-hyperelliptic
        FilterByALFixedPointsOnQuotient(~curves);

        // upward closure - if covering a non-hyperelliptic, then non-hyperelliptic
        UpwardClosure(~curves);

        // if a genus 3 covers a genus 2 curve, then it is hyperelliptic
        Genus3CoversGenus2(~curves);

        DownwardClosure(~curves);

        // Using Proposition 6 from [FH] adapted to the Shimura curve situation

        time FilterByComplicatedALFixedPointsOnQuotient(~curves); //long time

        UpwardClosure(~curves);

        // Using trace of Hecke operators to count points and show more curves are
        // non-hyperelliptic
        FilterByTrace(~curves);

        UpwardClosure(~curves);

        UpdateByIsomorphisms(~curves);

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

        FilterByWeilPolynomial(~curves);

        UpdateByIsomorphisms(~curves);

        UpwardClosure(~curves);

        DownwardClosure(~curves);

        Write("all_curves_progress.dat", "" : Overwrite);

        for curve in curves do
            Write("all_curves_progress.dat", Sprint(curve, "Magma"));
        end for;

    end if;

    if read_data then
        f :=  Read("all_curves_progress.dat");
        lines := Split(f, ";");
        curves := [ eval c : c in lines];
    end if;

    return curves;
end intrinsic;
