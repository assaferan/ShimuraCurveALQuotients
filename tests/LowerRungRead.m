// tests/LowerRungRead.m
//
// The solver script reads a two-cusp polytope rung off a higher cached rung with the same pole at 0:
// with only the enumerated (380, 93, 65) in a cache, a request for (380, 28, 65) is answered without
// Normaliz by the points of the higher rung whose pole at infinity is at most 28, and that set
// equals the enumerated (380, 28, 65).  A cached file carrying the derived-set mark is not used.
//
// The two rungs compared are the committed cache files, both Normaliz enumerations.

M := 380; nlo := 28; nhi := 93; m := 65;
lower := Sort(eval Read(Sprintf("polymake/polymake_solution_%o_%o_%o", M, nlo, m)));
tmp := Pipe("mktemp -d", "");
tmp := tmp[1..#tmp-1];
assert #tmp gt 0 and tmp[1] eq "/";
here := GetCurrentDirectory();
ok := System(Sprintf("mkdir -p %o/polymake && cp polymake/polymake_solution_%o_%o_%o %o/polymake/", tmp, M, nhi, m, tmp));
assert ok eq 0;
// The script says on its error stream which route answered a request: "lower rung:" for the
// read off a higher rung, a Normaliz line otherwise (Normaliz may or may not be installed here).
rc := System(Sprintf("cd %o && python3 %o/nmzsolve.py %o %o %o polymake/polymake_solution_%o_%o_%o 2> first.err",
                     tmp, here, M, nlo, m, M, nlo, m));
assert rc eq 0;
first_route := Read(Sprintf("%o/first.err", tmp));
got := Sort(eval Read(Sprintf("%o/polymake/polymake_solution_%o_%o_%o", tmp, M, nlo, m)));
// with the higher rung marked as derived, the read must not be used
ok := System(Sprintf("cd %o && rm polymake/polymake_solution_%o_%o_%o && (echo '// derived: test'; cat %o/polymake/polymake_solution_%o_%o_%o) > polymake/polymake_solution_%o_%o_%o",
                     tmp, M, nlo, m, here, M, nhi, m, M, nhi, m));
assert ok eq 0;
_ := System(Sprintf("cd %o && NMZ_TIMEOUT=60 python3 %o/nmzsolve.py %o %o %o polymake/polymake_solution_%o_%o_%o 2> second.err",
                    tmp, here, M, nlo, m, M, nlo, m));
second_route := Read(Sprintf("%o/second.err", tmp));
_ := System(Sprintf("rm -rf %o", tmp));

procedure test_LowerRungRead(lower, got, first_route, second_route)
    printf "Testing the lower-rung read, (380, 28, 65) off (380, 93, 65)...";
    assert "lower rung:" in first_route;                        // answered by the read, not by Normaliz
    assert got eq lower;
    assert not ("lower rung:" in second_route);                 // a derived higher rung is not read
    printf "Done! (%o points; a derived higher rung is not used)\n", #lower;
end procedure;

test_LowerRungRead(lower, got, first_route, second_route);
