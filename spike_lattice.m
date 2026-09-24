// Spike: inspect lattice structure of table targets + time the CM-point count.
// No Borcherds / EquationsOfCovers work -- just lattice queries + CM enumeration.
AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 0);

curves := GetHyperellipticCandidates();
printf "Loaded %o candidate curves.\n", #curves;

function neigh(curves, X)
    below := Sort([ <curves[i]`g, curves[i]`CurveID> : i in X`Covers ]);
    above := Sort([ <curves[i]`g, curves[i]`CurveID> : i in X`CoveredBy ]);
    return below, above;
end function;

procedure inspect(D, N, gens)
    printf "\n==== D=%o N=%o  target gens=%o ====\n", D, N, gens;
    if not exists(Xstar){X : X in curves | X`D eq D and X`N eq N and IsStarCurve(X)} then
        printf "  no star curve\n"; return;
    end if;
    b, a := neigh(curves, Xstar);
    printf "  Xstar CurveID=%o g=%o : Covers(below)=%o  CoveredBy(above, <g,id>)=%o\n",
        Xstar`CurveID, Xstar`g, b, a;
    W := AllALsFromGens(gens, D*N);
    if not exists(tj){j : j in [1..#curves] | curves[j]`D eq D and curves[j]`N eq N and curves[j]`W eq W} then
        printf "  TARGET W not found in candidate list at all\n"; return;
    end if;
    T := curves[tj];
    immediate := tj in Xstar`CoveredBy;
    bt, at := neigh(curves, T);
    printf "  TARGET CurveID=%o g=%o immediate-over-Xstar=%o\n", T`CurveID, T`g, immediate;
    printf "    Covers(below, <g,id>)=%o\n", bt;
    printf "    CoveredBy(above, <g,id>)=%o\n", at;
end procedure;

// Genus-2 Table 10 anchor (settled BIELLIPTIC)
inspect(6, 29, {3,29});
// Genus-0 Table 6 case known to work
inspect(10, 7, {2,5});
// Current failures (genus 2 targets)
inspect(34, 5, {2,5});
inspect(34, 5, {5,17});
inspect(26, 5, {2,5});

printf "\nSPIKE DONE.\n";
exit;
