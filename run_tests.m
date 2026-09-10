// usage: magma [target:=SUBSTRING] [exitsignal:=BOOL] [verbose:=INT] [debug:=BOOL] run_tests.m
if assigned filename then
  if "tests/" eq filename[1..6] then
    filename := filename[7..#filename];
  end if;
  tests := [filename];
else
  // ⚠ MIRROR THE CI MATRIX: .github/workflows/tests.yml uses grep -vE '^_|^run_filters\.m$'.
  // Files starting with "_" are shared helpers and hand-run tools, not tests; run_filters.m is the
  // full end-to-end pipeline (6h+, past GitHub's job limit) and is meant to be run deliberately.
  //
  // ⚠ WHY THIS MATTERS AND IS NOT COSMETIC: several helpers END WITH `exit;`, which terminates
  // Magma and SILENTLY TRUNCATES the suite. On 2026-09-09 `_gyinvol.m` did exactly that -- it
  // sorts AFTER the uppercase test names (ASCII "_" is 0x5F, above "Z"), so 34 X0_* tests ran,
  // then it called exit, and the remaining tests never ran and NO summary was printed. A truncated
  // run looks almost exactly like a clean one. The "_ is excluded" convention was documented in
  // tests/_crviso.m's header but had never actually been implemented here.
  tests := [f : f in Split(Pipe("cd tests && ls *.m && cd ..", ""), "\n")
            | f ne "" and f[1] ne "_" and f ne "run_filters.m"];
end if;
if assigned debug then
  SetDebugOnError(true);
end if;
AttachSpec("ShimuraQuotients.spec");

if assigned verbose then
  try
    verbose := StringToInteger(verbose);
  catch e
    verbose := 1;
  end try;
  SetVerbose("ShimuraQuotients", verbose);
end if;
failed := [];
if not assigned target then
  target := "";
end if;

counter := 0;
for filename in tests do
  if target in filename then
    counter +:=1;
    fullPath := "tests/" cat filename;
    timestamp := Time();
    if assigned debug then
      printf "%o: ", filename;
      assert eval (Read(fullPath) cat  "return true;");
      printf "Success! %o s\n", Time(timestamp);
    else
      try
        printf "%o: ", filename;
        assert eval (Read(fullPath) cat  "return true;");
        printf "Success! %o s\n", Time(timestamp);
      catch e
        Append(~failed, filename);
        printf "Fail! %o s\n %o\n", e, Time(timestamp);;
      end try;
    end if;
  end if;
end for;
if counter eq 0 then
  print "No matching target";
  exit 1;
end if;
if #failed gt 0 then
  print "Tests failed:";
  for f in failed do
    print f;
  end for;
end if;
if assigned exitsignal then
  exit #failed;
end if;

