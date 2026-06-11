// Fast class numbers of imaginary quadratic orders, backed by the precomputed
// LMFDB quadratic-imaginary class-group tables
// (https://www.lmfdb.org/NumberField/QuadraticImaginaryClassGroups).
//
// The tables only store *fundamental* discriminants.  For a fundamental
// discriminant d0 < 0 the class number h(d0) is looked up directly; for a
// non-fundamental discriminant D = d0 * f^2 the order class number is recovered
// from h(d0) via the conductor formula [Cox, Primes of the form x^2+ny^2, Thm 7.24]
//
//   h(D) = h(d0) * f / [O_K^* : O^*] * prod_{p | f} ( 1 - (d0/p)/p ).
//
// Data layout (see the LMFDB download page for the authoritative description):
//   * Negative fundamental discriminants are split by congruence of |d| into
//     four subdirectories:
//        cl3mod8  : |d| = 3 mod 8      cl7mod8  : |d| = 7 mod 8
//        cl4mod16 : |d| = 4 mod 16     cl8mod16 : |d| = 8 mod 16
//     Together these cover every negative fundamental discriminant.
//   * File k (cl{r}mod{m}.k.gz) covers k*2^28 <= |d| < (k+1)*2^28.
//   * Each line has three tab-separated fields  a  b  c1 c2 ... ct  where a is a
//     gap, b = h(d) is the class number, and the c_j are the invariant factors.
//     Reading lines in order, |d_i| = k*2^28 + r + m * (a_0 + a_1 + ... + a_i).
//
// Files are streamed once and cached on demand: a lookup for |d| only reads the
// table up to |d|, remembering every class number it passes, so subsequent
// lookups (which dominate) are served from memory.

CL_FILE_SPAN := 2^28;            // |d| range covered by a single data file
CL_DEFAULT_DIR := "/scratch/class-groups-quadratic-imaginary-fields";

// Module-level state, persistent across calls within a Magma session.
CL_STORE := NewStore();          // config + per-file stream cursors
CL_FUND := NewStore();           // |d0| -> h(d0)   for fundamental d0
CL_ORDER := NewStore();          // D    -> h(order of discriminant D)

function clGetAssoc(store)
    b, A := StoreIsDefined(store, "cache");
    if not b then A := AssociativeArray(); end if;
    return A;
end function;

function clDataDir()
    b, dir := StoreIsDefined(CL_STORE, "dir");
    if b then return dir; end if;
    env := GetEnv("CLASS_GROUPS_DIR");
    return (env ne "") select env else CL_DEFAULT_DIR;
end function;

intrinsic SetClassNumberDataDir(dir::MonStgElt)
{Set the directory that holds the precomputed quadratic-imaginary class-group
tables (overrides the CLASS_GROUPS_DIR environment variable and the default).}
    StoreSet(CL_STORE, "dir", dir);
end intrinsic;

intrinsic ClassNumberDataMaxAbsDisc() -> RngIntElt
{The exclusive upper bound on |d| covered by the precomputed class-group tables:
file indices run 0..4095 and file k covers k*2^28 <= |d| < (k+1)*2^28, so |d| < 2^40.
A discriminant of absolute value >= this bound is not in the tables and ClassNumberLU
falls back to a direct computation.}
    return 4096 * CL_FILE_SPAN;
end intrinsic;

// Congruence class <subdirectory/prefix, residue r, modulus m> of a negative
// fundamental discriminant with |d0| = m0.
function clResidue(m0)
    case m0 mod 8:
        when 3: return "cl3mod8", 3, 8;
        when 7: return "cl7mod8", 7, 8;
    end case;
    case m0 mod 16:
        when 4: return "cl4mod16", 4, 16;
        when 8: return "cl8mod16", 8, 16;
    end case;
    error Sprintf("|d| = %o is not a negative fundamental discriminant", m0);
end function;

// Look up h(d0) for a negative fundamental discriminant d0 with |d0| = m0,
// streaming the relevant table file on demand.  Returns false if the value is
// unavailable (missing data file / |d| beyond the downloaded range), in which
// case the caller should fall back to a direct computation.
function clFundClassNo(m0)
    fund := clGetAssoc(CL_FUND);
    b, h := IsDefined(fund, m0);
    if b then return true, h; end if;

    prefix, r, m := clResidue(m0);
    k := m0 div CL_FILE_SPAN;
    fkey := prefix cat "." cat IntegerToString(k);

    // Recover (or open) the streaming cursor for this file.
    have, st := StoreIsDefined(CL_STORE, fkey);
    if have and st`done then return false, _; end if;
    if not have then
        path := clDataDir() cat "/" cat prefix cat "/" cat fkey cat ".gz";
        if not FileExists(path) then
            // Mark exhausted so we don't probe again, and warn once.
            StoreSet(CL_STORE, fkey, rec<recformat<F, cum, done> | done := true>);
            wb, warned := StoreIsDefined(CL_STORE, "warned");
            if not (wb and warned) then
                printf "WARNING: class-group data file not found: %o\n", path;
                printf "         falling back to direct class number computation.\n";
                StoreSet(CL_STORE, "warned", true);
            end if;
            return false, _;
        end if;
        // 2>/dev/null suppresses zcat's harmless "stdout: Broken pipe" SIGPIPE
        // notice, which it prints whenever we stop reading the stream early (the
        // normal case: we break as soon as the target |d| is found, and abandon
        // the still-open pipe at process exit). The file's existence was already
        // checked above, so this does not mask a missing-file error.
        F := POpen("zcat " cat path cat " 2>/dev/null", "r");
        st := rec<recformat<F, cum, done> | F := F, cum := 0, done := false>;
    end if;

    base := k * CL_FILE_SPAN + r;   // |d| of the first line in this file
    cum := st`cum;
    F := st`F;
    found := false; result := 0;
    line := Gets(F);
    while not IsEof(line) do
        sp := Split(line, "\t");
        cum +:= StringToInteger(sp[1]);
        d := base + m * cum;
        hd := StringToInteger(sp[2]);
        fund[d] := hd;
        if d eq m0 then
            found := true; result := hd;
            line := "";                 // sentinel: stop, but cursor is mid-file
            break;
        elif d gt m0 then
            break;                       // passed it: m0 absent from this file
        end if;
        line := Gets(F);
    end while;

    st`cum := cum;
    if IsEof(line) then st`done := true; end if;
    StoreSet(CL_STORE, fkey, st);
    StoreSet(CL_FUND, "cache", fund);

    if found then return true, result; end if;
    return false, _;
end function;

intrinsic ClassNumberTableMaxDisc() -> RngIntElt
{Fundamental discriminants with |d0| at or above this value are served by a direct
ClassNumber computation instead of the streaming LMFDB tables.  The table lookup caches
every class number it streams past (~750 MB of memory per 2^28-sized data file, per
residue class), which is a large win for the small, frequently-reused discriminants that
dominate -- but pure cost for large discriminants, which are touched roughly once with no
reuse.  Capping the table range therefore bounds the per-process cache: at the default
2^28 it stays ~3 GB even when a high-genus point count reaches |d| ~ 10^10.  Tunable; raise
it to keep more lookups on the (faster, but memory-hungry) tables when memory permits.}
    return 2^28;
end intrinsic;

intrinsic ClassNumberLU(D::RngIntElt) -> RngIntElt
{The class number of the imaginary quadratic order of discriminant D (D < 0 and
D = 0 or 1 mod 4), using the precomputed LMFDB fundamental class-group tables for small
fundamental discriminants and a direct ClassNumber computation for large ones (see
ClassNumberTableMaxDisc) or whenever the tables do not cover D.}
    require D lt 0 and (D mod 4 in [0, 1]):
        "D must be a negative discriminant (D < 0 and D = 0 or 1 mod 4)";

    order := clGetAssoc(CL_ORDER);
    cached, h := IsDefined(order, D);
    if cached then return h; end if;

    D0 := FundamentalDiscriminant(D);

    // Large fundamental discriminants are served by a direct ClassNumber rather than the
    // streaming tables.  clFundClassNo caches every class number it streams past (~750 MB
    // per 2^28-sized file per residue class), so a single high-genus point count -- which
    // scatters lookups across |d| up to ~4*Qmax*p^g, i.e. tens of GB of files -- would
    // exhaust memory (this is what OOM'd FilterStarCurvesByFpAutomorphisms on a g=8 curve).
    // Such large discriminants are touched ~once with no reuse, so a direct computation is
    // both memory-safe and, in aggregate, cheaper than streaming+caching that whole range.
    if -D0 ge ClassNumberTableMaxDisc() then
        h := ClassNumber(D);
        order[D] := h;
        StoreSet(CL_ORDER, "cache", order);
        return h;
    end if;

    f := Isqrt(D div D0);                 // conductor; D = D0 * f^2

    ok, h0 := clFundClassNo(-D0);
    if not ok then
        h := ClassNumber(D);              // fallback (also covers non-maximal)
        order[D] := h;
        StoreSet(CL_ORDER, "cache", order);
        return h;
    end if;

    if f eq 1 then
        h := h0;
    else
        ui := 1;
        if   D0 eq -3 then ui := 3;
        elif D0 eq -4 then ui := 2;
        end if;
        prod := 1;
        for p in PrimeDivisors(f) do
            prod *:= (1 - KroneckerSymbol(D0, p) / p);
        end for;
        h := Integers()!(h0 * f / ui * prod);
    end if;

    order[D] := h;
    StoreSet(CL_ORDER, "cache", order);
    return h;
end intrinsic;

// Order class number h(D0*f^2) from the fundamental h(D0)  [Cox, Primes ..., Thm 7.24].
function clOrderFromFund(D0, f, h0)
    if f eq 1 then return h0; end if;
    ui := 1;
    if   D0 eq -3 then ui := 3;
    elif D0 eq -4 then ui := 2;
    end if;
    prod := 1;
    for p in PrimeDivisors(f) do
        prod *:= (1 - KroneckerSymbol(D0, p) / p);
    end for;
    return Integers()!(h0 * f / ui * prod);
end function;

// Batched fundamental class numbers.  m0s must be a sorted (ascending) sequence of positive
// |d0| for negative fundamental discriminants below the table range.  Streams each required
// data file exactly once, ascending, reading off only the requested discriminants -- so memory
// is O(#m0s), not O(#streamed).  Returns AssociativeArray |d0| :-> h(d0); any value the tables
// cannot supply (missing file, absent line) falls back to a direct ClassNumber.
function clFundClassNoBatch(m0s)
    res := AssociativeArray();
    byFile := AssociativeArray();              // <prefix, file k> :-> ascending [|d0|]
    for m0 in m0s do
        prefix, r, m := clResidue(m0);
        key := <prefix, m0 div CL_FILE_SPAN>;
        if not IsDefined(byFile, key) then byFile[key] := []; end if;
        Append(~byFile[key], m0);
    end for;
    for key in Keys(byFile) do
        prefix := key[1]; k := key[2];
        targets := byFile[key];
        path := clDataDir() cat "/" cat prefix cat "/" cat prefix cat "."
                cat IntegerToString(k) cat ".gz";
        if not FileExists(path) then
            for m0 in targets do res[m0] := ClassNumber(-m0); end for;
            continue;
        end if;
        _, r, m := clResidue(targets[1]);      // residue/modulus fixed within a prefix
        base := k * CL_FILE_SPAN + r;
        F := POpen("zcat " cat path cat " 2>/dev/null", "r");
        cum := 0; ti := 1;
        line := Gets(F);
        while (not IsEof(line)) and (ti le #targets) do
            sp := Split(line, "\t");
            cum +:= StringToInteger(sp[1]);
            d := base + m * cum;               // discriminants increase down the file
            while (ti le #targets) and (targets[ti] lt d) do
                res[targets[ti]] := ClassNumber(-targets[ti]);   // file skipped it: fall back
                ti +:= 1;
            end while;
            if (ti le #targets) and (d eq targets[ti]) then
                res[d] := StringToInteger(sp[2]);
                ti +:= 1;
            end if;
            line := Gets(F);
        end while;
        while ti le #targets do                // ran off the end of the file
            res[targets[ti]] := ClassNumber(-targets[ti]);
            ti +:= 1;
        end while;
    end for;
    return res;
end function;

intrinsic ClassNumberBatchLU(Ds::SeqEnum[RngIntElt]) -> Assoc
{Batched form of ClassNumberLU: returns an associative array D :-> h(D) for every D in Ds,
identical to calling ClassNumberLU(D) individually, but the table lookups are sorted by
fundamental discriminant and streamed in a single ascending pass per file, retaining only the
requested values.  Memory stays O(#Ds) and each file is read at most once, instead of caching
every class number streamed past -- use this for the whole discriminant set of a trace-formula
sum.  Large/uncovered fundamental discriminants use a direct ClassNumber, as in ClassNumberLU.}
    res := AssociativeArray();
    order := clGetAssoc(CL_ORDER);
    // The batch path streams request-only, so it is memory-safe over the whole table range
    // (unlike per-disc ClassNumberLU, which caches everything streamed and so caps out at
    // ClassNumberTableMaxDisc).  Only go direct for discriminants beyond the tables entirely.
    T := ClassNumberDataMaxAbsDisc();
    info := AssociativeArray();                // D :-> <D0, f>  for the table-served D
    needFund := {};                            // |d0| to fetch from the tables
    for D in Ds do
        require D lt 0 and (D mod 4 in [0,1]):
            "every D must be a negative discriminant (D < 0 and D = 0 or 1 mod 4)";
        if IsDefined(res, D) or IsDefined(info, D) then continue; end if;
        cached, h := IsDefined(order, D);
        if cached then res[D] := h; continue; end if;
        D0 := FundamentalDiscriminant(D);
        if -D0 ge T then
            h := ClassNumber(D);               // large fundamental: direct, memory-safe
            res[D] := h; order[D] := h;
        else
            info[D] := <D0, Isqrt(D div D0)>;
            Include(~needFund, -D0);
        end if;
    end for;
    h0map := clFundClassNoBatch(Sort(SetToSequence(needFund)));
    for D in Keys(info) do
        D0 := info[D][1]; f := info[D][2];
        h := clOrderFromFund(D0, f, h0map[-D0]);
        res[D] := h; order[D] := h;
    end for;
    StoreSet(CL_ORDER, "cache", order);
    return res;
end intrinsic;
