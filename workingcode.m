// This small procedure is used to update
// the curves after every filter by using
// upward and downward closures, as well as isomorphisms
intrinsic UpdateCurves(~curves::SeqEnum)
{Propagate hyperellipticity knowledge via isomorphisms and the coverage lattice.}
  UpdateByIsomorphisms(~curves);
  // upward closure - if X-->Y is dominant then gon(X) >= gon(Y).
  // [Poonen, A.1.(vii)]
  // In particular, if Y is non-hyperelliptic, so is X.
  UpwardClosure(~curves);
  // downward closure - if covered by subhyperelliptic, then subhyperelliptic
  DownwardClosure(~curves);
end intrinsic;

// The different stages of filtering
FILTER_STAGES := [*
    // <"FindPairs", FindPairs>,
    <"UpdateGenera", UpdateGenera>,
    // <"VerifyHHTable1", VerifyHHTable1>,
    <"UpdateByGenusStar", UpdateByGenus>,
    <"FilterByTraceStar", FilterByTrace>,
    // <"VerifyHHTable2", VerifyHHTable2>,
    <"HHProposition1", HHProposition1>,
    // <"VerifyHHProposition1", VerifyHHProposition1>,
    <"FilterStarCurvesByFpAutomorphisms", FilterStarCurvesByFpAutomorphisms>,
    // <"GetQuotientsAndGenera", GetQuotientsAndGenera>,
    <"UpdateByGenus", UpdateByGenus>,
    <"UpdateCurves1", UpdateCurves>,
    <"FilterByALFixedPointsOnQuotient", FilterByALFixedPointsOnQuotient>,
    <"UpdateCurves2", UpdateCurves>,
    <"Genus3CoversGenus2", Genus3CoversGenus2>,
    <"UpdateCurves3", UpdateCurves>,
    <"FilterByDegeneracyMorphism",FilterByDegeneracyMorphism>,
    <"UpdateCurves4", UpdateCurves>,
    <"FilterByComplicatedALFixedPointsOnQuotient", FilterByComplicatedALFixedPointsOnQuotient>,
    <"UpdateCurves5", UpdateCurves>,
    <"FilterByTrace", FilterByTrace>,
    <"UpdateCurves6", UpdateCurves>,
    <"FilterByWeilPolynomial",FilterByWeilPolynomialGenusScaled>,
    <"UpdateCurves7", UpdateCurves>,
    <"FilterByNonALInvolutions",FilterByNonALInvolutions>,
    <"UpdateCurves8", UpdateCurves>
*];

function compute_data(start_stage, stages)

    procedure run_stage(name, func, ~curves)
      t0 := Realtime();
      func(~curves);
      print name, " took ", Realtime() - t0;
      fname := Sprintf("data/curves_after_%o.dat", name);
      Write(fname, Sprint(curves, "Magma") : Overwrite);
      return;
    end procedure;

    // finding where to start
    start := 1;
    for i->stage in stages do
      if stage[1] eq start_stage then
        fname := Sprintf("data/curves_after_%o.dat", stage[1]);
        curves := eval Read(fname);
        start := i + 1;
      end if;
    end for;

    if (start eq 1) then
        // Find the largest prime we need to consider for the
        // inequality in Proposition 1.
        // TODO: Make this configurable
        r := GetLargestPrimeIndex();
        assert r eq 7;
        // Find all pairs (D,N) satisfying the inequality of
        // Proposition 1.
        t0 := Realtime();
        curves := FindPairs(r); // time (lava, 05/28/26) : 1.570
        print "FindPairs took ", Realtime() - t0;
        assert #curves eq 2342;
    end if;

    for i := start to #stages do
      stage := stages[i];
      // between FilterByFpAutomorphisms and UpdateByGenus,
      // we need to generate all curves from the star curves
      if (stage[1] eq "UpdateByGenus") then
	      t0 := Realtime();
        curves := GetQuotientsAndGenera(curves);
        // time (lava, 05/28/26) : 131.360
        print "GetQuotientsAndGenera took ", Realtime() - t0;
      end if;

      // run filtering stage
      run_stage(stage[1], stage[2], ~curves);

      // in certain cases, we add verifications
      case stage[1]:
        when "UpdateGenera":
	        VerifyHHTable1(curves);
        when "FilterByTraceStar":
	        VerifyHHTable2(curves);
        when "HHProposition1":
    	    VerifyHHProposition1(curves);
      end case;
   end for;

   return curves;
end function;

intrinsic GetHyperellipticCandidates(:recompute_data:=false,
				     read_data:=true,
				     start_stage:="",
             read_stage:="",
				     stages:=FILTER_STAGES) -> SeqEnum
{.}
    if read_data then
      if read_stage eq "" then read_stage := stages[#stages][1]; end if;
      fname := Sprintf("data/curves_after_%o.dat", read_stage);
      read_curves := eval Read(fname);
    end if;

    if recompute_data then
      curves := compute_data(start_stage, stages);
      assert #read_curves eq #curves;
      assert &and[(read_curves[j] eq curves[j]) : j in [1..#curves]];
    else
      curves := read_curves;
    end if;

    return curves;
end intrinsic;

// Last timing: 334554 s = 92.93 h = 3.87 d
