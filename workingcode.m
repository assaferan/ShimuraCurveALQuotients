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
        t0 := Realtime();
        star_curves := FindPairs(r); // time (lava, 05/28/26) : 1.570
        print "FindPairs took ", Realtime() - t0;
        // I added some code that just
        // focuses on the star quotients X_0^*(D,N)
        assert #star_curves eq 2342;
        t0 := Realtime();
        UpdateGenera(~star_curves); // time (lava, 05/28/26) : 0.900
        print "UpdateGenera for star curves took ", Realtime() - t0;
        VerifyHHTable1(star_curves);
        UpdateByGenus(~star_curves);
        t0 := Realtime();
        FilterByTrace(~star_curves); // time (lava, 05/28/26) : 6124.790
        print "FilterByTrace for star curves took ", Realtime() - t0;
        Write("data/star_curves_point_count.dat", Sprint(star_curves, "Magma") : Overwrite);
        VerifyHHTable2(star_curves);
        HHProposition1(~star_curves);
        VerifyHHProposition1(star_curves);
        t0 := Realtime();
        FilterStarCurvesByFpAutomorphisms(~star_curves, 10, 20); // time (lava, 05/28/26) : 8919.220
        print "FilterStarCurvesByFpAutomorphisms took ", Realtime() - t0;
        // Create a list of all Atkin-Lehner quotients
        // compute their genera, and store the covering structure.
        // writing to a file, in case we would like to load it directly
        Write("data/star_curves_after_FpAuto.dat", Sprint(star_curves, "Magma") : Overwrite);
        read_curves := eval Read("data/star_curves_after_FpAuto.dat");
        assert #read_curves eq #star_curves;
        assert &and[(read_curves[j] eq star_curves[j]) : j in [1..#star_curves]];
        t0 := Realtime();
        curves := GetQuotientsAndGenera(star_curves); // time (lava, 05/28/26) : 131.360
        print "GetQuotientsAndGenera took ", Realtime() - t0;
        // updating classification from the genera we computed
        t0 := Realtime();
        UpdateByGenus(~curves); // time (lava, 05/28/26) : 0.030
        print "UpdateByGenus took ", Realtime() - t0;
        UpdateCurves(~curves);
        Write("data/all_curves_after_genus.dat", Sprint(curves, "Magma") : Overwrite);
        // Using the fact that if w acts non-trivially and has more than
        // 4 fixed points on X, then X is non-hyperelliptic
        t0 := Realtime();
        FilterByALFixedPointsOnQuotient(~curves);
        print "FilterByALFixedPointsOnQuotient took ", Realtime() - t0;
        UpdateCurves(~curves);
        Write("data/all_curves_after_ALFixed.dat", Sprint(curves, "Magma") : Overwrite);
        // if a genus 3 covers a genus 2 curve, then it is hyperelliptic
        t0 := Realtime();
        Genus3CoversGenus2(~curves);
        print "Genus3CoversGenus2 took ", Realtime() - t0;
        UpdateCurves(~curves);
        Write("data/all_curves_after_genus3.dat", Sprint(curves, "Magma") : Overwrite);
        t0 := Realtime();
        FilterByDegeneracyMorphism(~curves);
        print "FilterByDegeneracyMorphism took ", Realtime() - t0;
        UpdateCurves(~curves);
        Write("data/all_curves_after_degeneracy.dat", Sprint(curves, "Magma") : Overwrite);
        // Using Proposition 6 from [FH] adapted to the Shimura curve situation
        t0 := Realtime();
        FilterByComplicatedALFixedPointsOnQuotient(~curves); // 64.840
        print "FilterByComplicatedALFixedPointsOnQuotient took ", Realtime() - t0;
        UpdateCurves(~curves);
        Write("data/all_curves_after_ComplicatedALFixed.dat", Sprint(curves, "Magma") : Overwrite);

        t0 := Realtime();
        FilterByNonALInvolutions(~curves);
        print "FilterByNonALInvolutions took ", Realtime() - t0;      
        UpdateCurves(~curves);
        Write("data/all_curves_after_NonALInvolutions.dat", Sprint(curves, "Magma") : Overwrite);

        // Using trace of Hecke operators to count points and show more curves are
        // non-hyperelliptic
        t0 := Realtime();
        FilterByTrace(~curves);
        print "FilterByTrace took ", Realtime() - t0;
        t0 := Realtime();
        UpdateCurves(~curves); // 13.750
        print "UpdateCurves took ", Realtime() - t0;
        Write("data/all_curves_after_trace.dat", Sprint(curves, "Magma") : Overwrite);
        t0 := Realtime();
        FilterByWeilPolynomial(~curves);
        print "FilterByWeilPolynomial took ", Realtime() - t0;
        UpdateCurves(~curves);
        Write("data/all_curves_after_Weil.dat", Sprint(curves, "Magma") : Overwrite);
    end if;

    if read_data then
        read_curves := eval Read("data/all_curves_after_Weil.dat");
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
