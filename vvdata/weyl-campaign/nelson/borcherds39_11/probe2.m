SetColumns(0);
AttachSpec("ShimuraQuotients.spec");
curves := GetHyperellipticCandidates();
// star genus for every base that has a committed model file
files := Split(Pipe("ls data/models/models_*.m", ""), "\n");
nzero := 0; npos := 0; posl := [];
for fn in files do
  s := Split(fn, "/");  s := s[#s];
  s := Substring(s, 8, #s-9); // strip models_ and .m
  parts := Split(s, "_");
  if #parts ne 2 then printf "SKIP %o\n", fn; continue; end if;
  D := StringToInteger(parts[1]); N := StringToInteger(parts[2]);
  ok := exists(Xs){X : X in curves | X`D eq D and X`N eq N and IsStarCurve(X)};
  if not ok then printf "NO STAR %o_%o\n", D,N; continue; end if;
  if Xs`g eq 0 then nzero +:= 1; else npos +:= 1; Append(~posl, <D,N,Xs`g>); end if;
end for;
printf "model bases with star genus 0 : %o\n", nzero;
printf "model bases with star genus >0: %o   %o\n", npos, posl;
// also: how many star curves overall have genus 0 vs >0
printf "DONE\n";
exit;
