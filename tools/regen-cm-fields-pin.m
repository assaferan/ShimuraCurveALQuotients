// tools/regen-cm-fields-pin.m -- regenerate the CM field-of-definition REGRESSION PIN
// (data/cm_fields_pin.m) from the CURRENT code only.
//
//   magma -b out:=/tmp/cm_fields_pin.m tools/regen-cm-fields-pin.m < /dev/null     # from repo root
//   magma -b out:=/tmp/x.m pin:=data/cm_fields_pin.m tools/regen-cm-fields-pin.m < /dev/null
//
// `out:=` is REQUIRED and there is deliberately no default: this script must never overwrite the
// committed pin by accident.  `pin:=` (default data/cm_fields_pin.m) is the committed file whose
// key set -- PIN_CURVES, PIN_DISCS and PIN_EXCLUDED -- is reused unchanged, so the regenerated file
// has exactly the same keys in the same order and a plain diff shows only changed VALUES.
//
// THEN DIFF AGAINST THE COMMITTED FILE, IGNORING THE HEADER COMMENT:
//   diff <(grep -v '^//' data/cm_fields_pin.m) <(grep -v '^//' /tmp/cm_fields_pin.m)
// and REVIEW EVERY CHANGED LINE.  A changed entry means the code now gives a different field of
// definition (or a different existence verdict) for that CM point than it did when pinned.  That
// is either a bug you just introduced or a fix you intend; this script cannot tell which.
// Copying the output over data/cm_fields_pin.m ("re-pinning") is a statement that you INTEND every
// such change -- say why in the commit message, and check changed points against a source
// (tests/CMPoints.m, tests/CMFieldsOfDefinition.m) where one exists.
//
// NOTE ON PROVENANCE.  The committed pin was made from THREE versions of the code (081d1ae,
// fddce4e, working tree of 2026-09-24) and keeps only values they all agree on.  This script sees
// only the current code, so its output is a weaker object; its header says so.  If you re-pin,
// keep that distinction honest in the header.
//
// Per key it computes FieldsOfDefinitionOfCMPoint (slow), FieldsOfDefinitionOfCMPointFast and
// DegreeOfFieldOfDefinitionOfCMPoint.  A key where slow and fast differ up to isomorphism, either
// errors, or the degree function disagrees with the field degrees is reported as INCONSISTENT and
// OMITTED from PIN (so it shows up in the diff as a deleted line).  Runtime ~6 min.
if not assigned out then
    print "usage: magma -b out:=OUTFILE [pin:=data/cm_fields_pin.m] tools/regen-cm-fields-pin.m < /dev/null";
    print "ERROR: out:= is required (no default, so the committed pin is never overwritten by accident).";
    exit 1;
end if;
if not assigned pin then pin := "data/cm_fields_pin.m"; end if;
SetColumns(0);   // no line wrapping in the written file; do not rely on a local .magmarc
AttachSpec("ShimuraQuotients.spec");

src := Read(pin);
PIN_CURVES   := eval (src cat "\nreturn PIN_CURVES;");
PIN_DISCS    := eval (src cat "\nreturn PIN_DISCS;");
PIN_EXCLUDED := eval (src cat "\nreturn PIN_EXCLUDED;");
excluded := {<e[1], e[2], e[3], e[4]> : e in PIN_EXCLUDED};
PIN          := eval (src cat "\nreturn PIN;");
committed := AssociativeArray();
for e in PIN do committed[<e[1], e[2], e[3], e[4]>] := e[5]; end for;

// Description of a list of fields as a sorted set of coefficient lists of absolute defining
// polynomials; [0, 1] is Q.  Built-in Magma only -- CI's Magma has no Polredabs (that comes from
// a locally attached package).  So that a diff shows only REAL changes, a field isomorphic to one
// of the committed entry's fields (refs) is written with the committed polynomial; a genuinely new
// field is written as DefiningPolynomial(OptimizedRepresentation(.)).
nf := func<c | c eq [0, 1] select Rationals() else NumberField(Polynomial(Rationals(), c))>;
iso := func<F, G | AbsoluteDegree(F) eq AbsoluteDegree(G) and
                   (AbsoluteDegree(F) eq 1 or IsIsomorphic(AbsoluteField(F), AbsoluteField(G)))>;
canon := function(Fs, refs)
    S := {};
    for F in Fs do
        if Type(F) eq FldRat or AbsoluteDegree(F) eq 1 then
            Include(~S, [0, 1]); continue;
        end if;
        if exists(c){c : c in refs | iso(nf(c), F)} then
            Include(~S, c);
        else
            Include(~S, Coefficients(DefiningPolynomial(OptimizedRepresentation(AbsoluteField(F)))));
        end if;
    end for;
    return Sort(Setseq(S));
end function;
// equal as sets up to isomorphism
same_fields := function(A, B)
    return &and[exists{G : G in B | iso(F, G)} : F in A] and &and[exists{F : F in A | iso(F, G)} : G in B];
end function;
// Format exactly as data/cm_fields_pin.m (no Magma line wrapping).
seqstr := func<s | "[" cat Join([Sprint(x) : x in s], ", ") cat "]">;
fldstr := func<L | "[" cat Join([seqstr(c) : c in L], ", ") cat "]">;

entries := [];
bad := 0;
t0 := Cputime();
for c in PIN_CURVES do
    D, N, W := Explode(c);
    X := CreateShimuraQuot(D, N, Set(W));
    for d in PIN_DISCS do
        if <D, N, W, d> in excluded then continue; end if;
        key := Sprintf("%o, %o, %o, %o", D, N, seqstr(W), d);
        refs := IsDefined(committed, <D, N, W, d>) select committed[<D, N, W, d>] else [];
        try
            slowF := FieldsOfDefinitionOfCMPoint(X, d);
        catch e
            printf "INCONSISTENT <%o>: slow raised %o\n", key, e`Object; bad +:= 1; continue;
        end try;
        try
            fastF := FieldsOfDefinitionOfCMPointFast(X, d);
        catch e
            printf "INCONSISTENT <%o>: fast raised %o\n", key, e`Object; bad +:= 1; continue;
        end try;
        deg := DegreeOfFieldOfDefinitionOfCMPoint(X, d);
        fast := canon(fastF, refs);
        degs := {#p - 1 : p in fast};
        if not same_fields(slowF, fastF) then
            printf "INCONSISTENT <%o>: slow %o fast %o\n", key, canon(slowF, refs), fast; bad +:= 1; continue;
        end if;
        if degs ne (deg eq 0 select {} else {deg}) then
            printf "INCONSISTENT <%o>: degree function %o, field degrees %o\n", key, deg, degs;
            bad +:= 1; continue;
        end if;
        Append(~entries, "<" cat key cat ", " cat fldstr(fast) cat ", " cat Sprint(deg) cat ">");
    end for;
    printf "%o %o %o done (%o s so far)\n", D, N, W, Cputime(t0);
end for;

today := Pipe("date +%Y-%m-%d", "");
today := today[1..#today-1];
header := "// data/cm_fields_pin.m -- REGRESSION PIN of fields of definition of CM points.\n"
  cat "//\n"
  cat "// REGENERATED by tools/regen-cm-fields-pin.m from the CURRENT code only, on " cat today cat ",\n"
  cat "// reusing the key set of " cat pin cat ".  This is the code's output, NOT a source of truth,\n"
  cat "// and unlike the committed pin it is NOT cross-checked against older versions.\n"
  cat "// Diff it against the committed file and review every change before re-pinning.\n\n";

curvestr := Join([Sprintf("<%o, %o, %o>", c[1], c[2], seqstr(c[3])) : c in PIN_CURVES], ",\n  ");
exclstr := Join([Sprintf("<%o, %o, %o, %o, \"%o\">", e[1], e[2], seqstr(e[3]), e[4], e[5])
                 : e in PIN_EXCLUDED], ",\n");
body := "PIN_CURVES := [\n  " cat curvestr cat "\n];\n\nPIN_DISCS := " cat seqstr(PIN_DISCS) cat ";\n\n"
  cat "PIN_EXCLUDED := [\n" cat exclstr cat "\n];\n\n"
  cat "PIN := [*\n" cat Join(entries, ",\n") cat "\n*];\n";
Write(out, header cat body : Overwrite := true);
printf "wrote %o: %o entries, %o INCONSISTENT keys omitted, %o excluded keys skipped (%o s)\n",
    out, #entries, bad, #excluded, Cputime(t0);
quit;
