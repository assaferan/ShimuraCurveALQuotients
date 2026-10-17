SetColumns(0);
AttachSpec("ShimuraQuotients.spec");
curves := GetHyperellipticCandidates();
stars := [X : X in curves | IsStarCurve(X) and X`D ne 1];
printf "#star curves (D>1) = %o\n", #stars;
g0 := [X : X in stars | X`g eq 0];
printf "  genus 0: %o ; genus >=1: %o\n", #g0, #stars-#g0;
// subhyperelliptic quotients sitting over a positive-genus star
nsub_pos := 0; nsub_zero := 0;
for X in stars do
  for i in X`CoveredBy do
    Y := curves[i];
    if assigned Y`IsSubhyp and Y`IsSubhyp then
      if X`g eq 0 then nsub_zero +:= 1; else nsub_pos +:= 1; end if;
    end if;
  end for;
end for;
printf "subhyperelliptic immediate covers: over genus-0 star %o ; over positive-genus star %o\n", nsub_zero, nsub_pos;
// 39_11 target
assert exists(T){X : X in curves | X`CurveID eq 9679};
printf "target 9679: IsSubhyp assigned=%o", assigned T`IsSubhyp;
if assigned T`IsSubhyp then printf " value=%o", T`IsSubhyp; end if;
printf "  IsHyp assigned=%o", assigned T`IsHyp;
if assigned T`IsHyp then printf " value=%o", T`IsHyp; end if;
printf "\n  TestInWhichProved: %o\n", assigned T`TestInWhichProved select T`TestInWhichProved else "(none)";
printf "DONE\n";
exit;
