// tests/_modelfile.m -- read a committed model set, isolated from any enclosing scope.
//
// ⚠ WHY THIS IS ITS OWN FILE. Reading data/models/models_D_N.m needs `eval` (Magma's `load` takes
// a literal filename and so cannot iterate or take D, N as arguments), and an `eval` inside a
// procedure that ALSO CLOSES OVER AN OUTER VARIABLE segfaults Magma 2.29 -- the same trap that
// forces ModelChecks.m and ModelRegen.m to be written as top-level statements. This function
// closes over nothing, so the eval is safe here and callers get a plain associative array back.
//
// Returns: found (BoolElt), models (Assoc keyed by Sort(W) as a sequence of integers).

function ReadModelSet(D, N)
    // ⚠ NOT `Pipe("ls ...")` to test for the file: Pipe raises on a non-zero exit status, so the
    // absent-file case threw instead of returning false. Read's own failure is the test.
    fname := Sprintf("data/models/models_%o_%o.m", D, N);
    ok := true; txt := "";
    try txt := Read(fname); catch e ok := false; end try;
    if not ok then return false, AssociativeArray(); end if;
    return true, eval (txt cat "\nreturn models;");
end function;
