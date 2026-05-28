// This small procedure is used to update
// the curves after every filter by using
// upward and downward closures, as well as isomorphisms
procedure UpdateCurves(~curves)
  UpdateByIsomorphisms(~curves);
  // upward closure - if X-->Y is dominant then gon(X) >= gon(Y).
  // [Poonen, A.1.(vii)]
  // In particular, if Y is non-hyperelliptic, so is X.
  UpwardClosure(~curves);
  // downward closure - if covered by subhyperelliptic, then subhyperelliptic
  DownwardClosure(~curves);
end procedure;

intrinsic GetHyperellipticCandidates(:recompute_data:=false, read_data :=true) -> SeqEnum
{.}
    // SetDebugOnError(true);
    // SetVerbose("ShimuraQuotients", 3);

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
        FilterByTrace(~star_curves); // time : 2850.270
        VerifyHHTable2(star_curves);
        // Create a list of all Atkin-Lehner quotients
        // compute their genera, and store the covering structure.
        // writing to a file, in case we would like to load it directly
        Write("data/star_curves_point_count.dat", Sprint(star_curves, "Magma") : Overwrite);
        read_curves := eval Read("data/star_curves_point_count.dat");
        assert #read_curves eq #star_curves;
        assert &and[(read_curves[j] eq star_curves[j]) : j in [1..#star_curves]];
        HHProposition1(~star_curves);
        VerifyHHProposition1(star_curves);
        FilterStarCurvesByFpAutomorphisms(~star_curves, 10, 20);

        time curves := GetQuotientsAndGenera(star_curves); // 61.520
        // updating classification from the genera we computed
        UpdateByGenus(~curves);
        UpdateCurves(~curves);
        // Using the fact that if w acts non-trivially and has more than
        // 4 fixed points on X, then X is non-hyperelliptic
        FilterByALFixedPointsOnQuotient(~curves);
        UpdateCurves(~curves);
        // if a genus 3 covers a genus 2 curve, then it is hyperelliptic
        Genus3CoversGenus2(~curves);
        UpdateCurves(~curves);
        // Using Proposition 6 from [FH] adapted to the Shimura curve situation
        time FilterByComplicatedALFixedPointsOnQuotient(~curves); // 64.840
        UpdateCurves(~curves);
        // Using trace of Hecke operators to count points and show more curves are
        // non-hyperelliptic
        time FilterByTrace(~curves);
        time UpdateCurves(~curves); // 13.750
        Write("data/all_curves_progress.dat", Sprint(curves, "Magma") : Overwrite);
        Write("data/all_curves_progress.dat", Sprint(curves, "Magma") : Overwrite);
        UpdateCurves(~curves);
        Write("data/all_curves_progress.dat", Sprint(curves, "Magma") : Overwrite);
        time FilterByWeilPolynomial(~curves);
        UpdateCurves(~curves);
        Write("data/all_curves_progress.dat", Sprint(curves, "Magma") : Overwrite);
        time FilterByDegeneracyMorphism(~curves);
        UpdateCurves(~curves);
        Write("data/all_curves_progress.dat", Sprint(curves, "Magma") : Overwrite);
    end if;

    if read_data then
        read_curves := eval Read("data/all_curves_progress.dat");
        if recompute_data then
            assert #read_curves eq #curves;
            assert &and[(read_curves[j] eq curves[j]) : j in [1..#curves]];
        else
            curves := read_curves;
        end if;
    end if;

    return curves;
end intrinsic;

// Last timing: 334554 s = 92.93 h = 3.87 d
