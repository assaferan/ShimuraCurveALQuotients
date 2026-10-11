AttachSpec("ShimuraQuotients.spec");
Attach("/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/fod/FodLambda.m");
SetColumns(0);
load "/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/fod/gydata.m";
curves := GetHyperellipticCandidates();
nchecked := 0; nambig := 0; nambig_ok := 0; nmismatch := 0; nunique_disagree := 0;
for datum in data do
    D, N, W, cm := Explode(datum);
    X := rep{X : X in curves | X`D eq D and X`N eq N and X`W eq W};
    for pt in cm do
        d, K := Explode(pt);
        if ClassNumber(d) gt 8 then continue; end if;
        if AbsoluteDegree(K) eq 1 then continue; end if;      // degree 1: nothing to choose
        Fs, nvalid := FieldsOfDefinitionOfCMPointLambda(X, d);
        old := FieldsOfDefinitionOfCMPointFast(X, d);
        if #Fs eq 0 then printf "NO POINT (D,N,W,d)=(%o,%o,%o,%o)\n", D, N, W, d; continue; end if;
        nchecked +:= 1;
        okK := IsIsomorphic(Fs[1], K);
        disputed := exists{g : g in degree_only_pts | g eq <D, N, W, d>};
        if nvalid ge 2 then
            nambig +:= 1;
            if okK then nambig_ok +:= 1; end if;
            printf "AMBIGUOUS (%o,%o,%o,d=%o): Fast gave %o fields, Lambda picks %o, GY field %o -> %o%o\n", D, N, W, d, #old, DefiningPolynomial(Fs[1]), DefiningPolynomial(K), okK select "MATCH" else "MISMATCH", disputed select " (disputed in test)" else "";
        else
            if not okK then nmismatch +:= 1; printf "UNIQUE MISMATCH (%o,%o,%o,d=%o): Lambda %o vs GY %o%o\n", D, N, W, d, DefiningPolynomial(Fs[1]), DefiningPolynomial(K), disputed select " (disputed in test)" else ""; end if;
            if not (#old eq 1 and IsIsomorphic(old[1], Fs[1])) then nunique_disagree +:= 1; printf "UNIQUE but Lambda != Fast (%o,%o,%o,d=%o)\n", D, N, W, d; end if;
        end if;
    end for;
end for;
printf "SUMMARY: %o degree>=2 GY points checked; algebra test ambiguous at %o of them, Lambda matches GY at %o of those; unique cases mismatching GY: %o; unique cases where Lambda != Fast: %o\n", nchecked, nambig, nambig_ok, nmismatch, nunique_disagree;
