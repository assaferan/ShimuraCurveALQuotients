// tests/NormalizOrder.m
//
// The solver script hands Normaliz its coordinates in one of two orders (the divisors of M, or the
// divisors farthest from sqrt(M) first) and writes the solution back in divisor order.  Both orders
// must give the same set of lattice points, and the same set as the committed cache file.  The
// polytope is (60, 11, 0), 163 points, a few seconds of Normaliz; NORMALIZ_BIN must point to a
// Normaliz binary (CI has /usr/bin/normaliz).
//
// The reference is the committed cache file polymake/polymake_solution_60_11_0, a Normaliz
// enumeration from before the order was configurable.

M := 60; n := 11;
cached := Sort(eval Read(Sprintf("polymake/polymake_solution_%o_%o_0", M, n)));
here := GetCurrentDirectory();
got := AssociativeArray();
for order in ["ends", "id"] do
    tmp := Pipe("mktemp -d", "");
    tmp := tmp[1..#tmp-1];
    assert #tmp gt 0 and tmp[1] eq "/";
    ok := System(Sprintf("mkdir -p %o/polymake", tmp));          // an empty cache: the solve goes to Normaliz
    assert ok eq 0;
    rc := System(Sprintf("cd %o && NMZ_ORDER=%o python3 %o/nmzsolve.py %o %o 0 polymake/polymake_solution_%o_%o_0 2> nmzsolve.err",
                         tmp, order, here, M, n, M, n));
    assert rc eq 0;
    got[order] := Sort(eval Read(Sprintf("%o/polymake/polymake_solution_%o_%o_0", tmp, M, n)));
    _ := System(Sprintf("rm -rf %o", tmp));
end for;

procedure test_NormalizOrder(cached, got)
    printf "Testing that both coordinate orders give the lattice points of (60, 11, 0)...";
    assert got["ends"] eq cached;
    assert got["id"] eq cached;
    printf "Done! (%o points in both orders)\n", #cached;
end procedure;

test_NormalizOrder(cached, got);
