// Model generation with a GENUS CAP on the targets, plus IntegralSolution.
//
//   NORMALIZ_BIN=... magma -b D_s:=74 N_s:=3 GCAP:=2 OUTDIR:=/path gencapped.m
//
// CM-point demand is num_vals = max(2g+5) over the RETAINED covers, so one high-genus sibling
// inflates the requirement for all of them.  cmsupply.m says 74_3, 74_5 and 58_5 are SHORT only
// because of their top-genus cover: capping the genus drops the demand below the available
// supply while still keeping most covers.  This restricts Targets to the covers with g <= GCAP.
SetColumns(0);
// Progress is ON by default.  A silent multi-hour run is worthless when the question is
// "is it progressing or stuck" -- the [n/6] stage markers in EquationsCovers.m are vprintf
// at level 1, so without this the log stays empty until the run ends.  Set VERB:=0 to mute.
AttachSpec("ShimuraQuotients.spec");
verb := 1;
if assigned VERB then verb := StringToInteger(VERB); end if;
SetVerbose("ShimuraQuotients", verb);
Px<x> := PolynomialRing(Rationals());
D := StringToInteger(D_s); N := StringToInteger(N_s);
gcap := StringToInteger(GCAP);
outdir := "data/models";
if assigned OUTDIR then outdir := OUTDIR; end if;

curves := GetHyperellipticCandidates();
Xstar := rep{X : X in curves | X`D eq D and X`N eq N and IsStarCurve(X)};

all_g := [curves[i]`g : i in Xstar`CoveredBy];
keep  := [i : i in Xstar`CoveredBy | curves[i]`g le gcap];
Targets := {curves[i]`W : i in keep};
printf "GENCAP base %o_%o gcap %o : genera %o -> keeping %o of %o covers, demand %o -> %o\n",
       D, N, gcap, Sort(all_g), #keep, #Xstar`CoveredBy,
       Maximum([2*g+5 : g in all_g]),
       IsEmpty(keep) select 0 else Maximum([2*curves[i]`g+5 : i in keep]);
if IsEmpty(keep) then printf "GENCAP RESULT %o_%o EMPTY\n", D, N; quit; end if;

t0 := Realtime();
try
    covers, ws := AllEquationsAboveCovers(Xstar, curves : Prec := 100,
                                          IntegralSolution := true, Targets := Targets);
    printf "GENCAP RESULT %o_%o OK covers %o time %os\n", D, N, #Keys(covers), Realtime()-t0;
    agg := AssociativeArray();
    for label in Keys(covers) do
      X := curves[label]; Wkey := Sort([w : w in X`W]);
      if not IsDefined(agg, Wkey) then agg[Wkey] := [* *]; end if;
      for base in Keys(covers[label]) do Append(~agg[Wkey], <X`g, covers[label][base]>); end for;
    end for;
    lines := [ Sprintf("// Genus-capped (g <= %o) cover models for X_0(%o,%o)*, IntegralSolution", gcap, D, N),
               "P<x> := PolynomialRing(Rationals());", "models := AssociativeArray();" ];
    for Wkey in Keys(agg) do
      for e in agg[Wkey] do
        Append(~lines, Sprintf("models[%o] := [* <%o, %o> *];", Wkey, e[1], e[2]));
      end for;
    end for;
    Write(Sprintf("%o/capped_%o_%o.m", outdir, D, N), Join(lines, "\n") : Overwrite);
    printf "GENCAP WROTE %o/capped_%o_%o.m (%o keys)\n", outdir, D, N, #Keys(agg);
catch e
    obj := Sprint(e`Object);
    printf "GENCAP RESULT %o_%o FAIL time %os :: %o\n", D, N, Realtime()-t0,
           Join(Split((#obj gt 120) select obj[1..120] else obj, "\n"), " | ");
end try;
printf "GENCAP DONE %o_%o\n", D, N;
quit;
