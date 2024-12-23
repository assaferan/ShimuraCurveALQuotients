AttachSpec("shimuraquots.spec");
curves := GetHyperellipticCandidates();
unknown := [c : c in curves | not assigned c`IsSubhyp];
unknownstar := [ c : c in unknown | IsStarCurve(c)];
unknownmodstar := [ c : c in unknownstar | c`D eq 1];
lut := AssociativeArray();
for i in [1..#curves] do
    c := curves[i];
    lut[<c`D, c`N, c`W>] := i;
end for;
unknowng4 := [c : c in curves | not assigned c`IsSubhyp and c`g eq 4];