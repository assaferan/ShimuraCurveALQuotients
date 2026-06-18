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

// ===========================================================================
// Fast random-access lookup.
//
// The streaming path above caches every class number it reads past, which makes
// repeated small-|d| lookups very cheap but costs ~750 MB of memory per 2^28-sized
// file -- so ClassNumberTableMaxDisc has to cap the table range to bound memory.
//
// The path below instead extracts each data file once into a fixed-width, sorted
// (|d|, h) record file and answers each lookup by binary search: O(log n) seeks and
// O(1) memory, no cache-everything.  That removes the memory pressure, so the tables
// can serve the entire downloaded range (|d| < 2^40) and the Weil disc-depth budget
// can be raised freely.  Extraction is one-time (~30 s/file) and persisted on disk
// (CLASS_GROUPS_FAST_DIR), shared across runs.
// ===========================================================================
CL_FAST_DW := 13;                            // |d| field width (|d| < 10^13 > 2^40)
CL_FAST_HW := 8;                             // h field width   (h   < 10^8)
CL_FAST_RW := CL_FAST_DW + CL_FAST_HW + 1;   // record width including the newline

function clFastDir()
    b, dir := StoreIsDefined(CL_STORE, "fastdir");
    if b then return dir; end if;
    env := GetEnv("CLASS_GROUPS_FAST_DIR");
    return (env ne "") select env else "/scratch/class-groups-fast";
end function;

intrinsic SetClassNumberFastDir(dir::MonStgElt)
{Set the directory holding the extracted fixed-width (|d|,h) record files used by the
fast binary-search class-number lookup (overrides the CLASS_GROUPS_FAST_DIR environment
variable and the built-in default).}
    StoreSet(CL_STORE, "fastdir", dir);
end intrinsic;

// Extract the fixed-width record file for (prefix, file k) from the gzipped LMFDB table,
// reconstructing |d| = base + m*(running gap sum) and writing "%013d%08d\n" per line (the
// source is already ascending in |d|).  Atomic via a temp file + rename, so concurrent
// workers are safe.  Returns the path on success, or "" if the source data file is missing.
function clExtractFastFile(prefix, k, r, m, fwdir, fwpath)
    gzpath := clDataDir() cat "/" cat prefix cat "/" cat prefix cat "."
              cat IntegerToString(k) cat ".gz";
    if not FileExists(gzpath) then return ""; end if;
    base := k * CL_FILE_SPAN + r;
    cmd := "mkdir -p " cat fwdir cat " && tmp=$(mktemp " cat fwdir cat "/.tmpXXXXXX) && zcat "
        cat gzpath cat " 2>/dev/null | awk -F'\\t' -v base=" cat IntegerToString(base)
        cat " -v m=" cat IntegerToString(m)
        cat " 'BEGIN{cum=0} {cum+=$1; d=base+m*cum; if($2>=100000000){print \"HOVERFLOW\">\"/dev/stderr\"; exit 1} printf \"%013d%08d\\n\", d, $2}' > $tmp && mv -f $tmp "
        cat fwpath;
    _ := System(cmd);
    if FileExists(fwpath) then return fwpath; end if;
    return "";
end function;

// Return <open handle, record count> for the fixed-width file of (prefix, file k), extracting
// it on first use.  Caches the handle + count in CL_STORE.  Returns false if data is missing.
function clFastFile(prefix, k, r, m)
    fkey := "fast:" cat prefix cat "." cat IntegerToString(k);
    have, st := StoreIsDefined(CL_STORE, fkey);
    if have then
        if st`missing then return false, _, _; end if;
        return true, st`F, st`nrec;
    end if;
    fwdir := clFastDir() cat "/" cat prefix;
    fwpath := fwdir cat "/" cat prefix cat "." cat IntegerToString(k) cat ".fw";
    if not FileExists(fwpath) then
        _ := clExtractFastFile(prefix, k, r, m, fwdir, fwpath);
    end if;
    if not FileExists(fwpath) then
        StoreSet(CL_STORE, fkey, rec<recformat<F, nrec, missing> | missing := true>);
        return false, _, _;
    end if;
    nrec := StringToInteger(Pipe("stat -c%s " cat fwpath cat " | tr -dc 0-9", "")) div CL_FAST_RW;
    F := Open(fwpath, "r");
    StoreSet(CL_STORE, fkey, rec<recformat<F, nrec, missing> | F := F, nrec := nrec, missing := false>);
    return true, F, nrec;
end function;

// h(d0) for a negative fundamental d0 with |d0| = m0, by binary search in the extracted file.
// O(log n) seeks, O(1) memory, no caching.  Returns false if unavailable (missing file / |d|
// beyond the downloaded range / line absent), so the caller can fall back to a direct compute.
function clFundClassNoFast(m0)
    prefix, r, m := clResidue(m0);
    k := m0 div CL_FILE_SPAN;
    ok, F, nrec := clFastFile(prefix, k, r, m);
    if not ok then return false, _; end if;
    lo := 0; hi := nrec - 1;
    while lo le hi do
        mid := (lo + hi) div 2;
        Seek(F, mid * CL_FAST_RW, 0);
        s := Read(F, CL_FAST_DW + CL_FAST_HW);
        d := StringToInteger(Substring(s, 1, CL_FAST_DW));
        if d eq m0 then
            return true, StringToInteger(Substring(s, CL_FAST_DW + 1, CL_FAST_HW));
        elif d lt m0 then lo := mid + 1;
        else hi := mid - 1;
        end if;
    end while;
    return false, _;
end function;

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
D = 0 or 1 mod 4).  Fundamental discriminants within the downloaded LMFDB table range
(|d0| < 2^40) are served by a binary-search lookup in the extracted tables (O(1) memory,
no cache-everything); larger ones, or any the tables do not cover, use a direct ClassNumber
computation.  The order class number is recovered from the fundamental one via the conductor
formula [Cox, Primes ..., Thm 7.24].}
    require D lt 0 and (D mod 4 in [0, 1]):
        "D must be a negative discriminant (D < 0 and D = 0 or 1 mod 4)";

    order := clGetAssoc(CL_ORDER);
    cached, h := IsDefined(order, D);
    if cached then return h; end if;

    D0 := FundamentalDiscriminant(D);

    // The fast lookup binary-searches an extracted file -- O(1) memory, no cache-everything --
    // so it serves the whole downloaded range without the memory blow-up that forced the
    // streaming path's ClassNumberTableMaxDisc cap (this is what OOM'd
    // FilterStarCurvesByFpAutomorphisms on a g=8 curve).  Only fall back to a direct
    // ClassNumber for fundamental discriminants beyond the tables entirely, or missing data.
    if -D0 ge ClassNumberDataMaxAbsDisc() then
        h := ClassNumber(D);
        order[D] := h;
        StoreSet(CL_ORDER, "cache", order);
        return h;
    end if;

    ok, h0 := clFundClassNoFast(-D0);
    if not ok then
        h := ClassNumber(D);              // fallback (missing data file / line absent)
        order[D] := h;
        StoreSet(CL_ORDER, "cache", order);
        return h;
    end if;

    f := Isqrt(D div D0);                 // conductor; D = D0 * f^2
    h := clOrderFromFund(D0, f, h0);
    order[D] := h;
    StoreSet(CL_ORDER, "cache", order);
    return h;
end intrinsic;

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
