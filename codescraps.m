AttachSpec("ShimuraQuotients.spec");
curves := GetHyperellipticCandidates();
unknown := [c : c in curves | not assigned c`IsSubhyp];
unknownmod := [c : c in unknown | c`D eq 1];
lut := AssociativeArray();
for i in [1..#curves] do
    c := curves[i];
    lut[<c`D, c`N, c`W>] := i;
end for;
unknownstar := [ c : c in unknown | IsStarCurve(c)];
unknownmodstar := [ c : c in unknownstar | c`D eq 1];
unknowng3 := [c : c in curves | not assigned c`IsSubhyp and c`g eq 3];
unknowng4 := [c : c in curves | not assigned c`IsSubhyp and c`g eq 4];
{* u`g : u in unknown *};
unknowng5 := [c : c in curves | not assigned c`IsSubhyp and c`g eq 5];
genera := {*  c`g : c in unknown*};
generastar := {*  c`g : c in unknownstar*};