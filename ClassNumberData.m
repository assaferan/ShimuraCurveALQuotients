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
        F := POpen("zcat " cat path, "r");
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

intrinsic ClassNumberLU(D::RngIntElt) -> RngIntElt
{The class number of the imaginary quadratic order of discriminant D (D < 0 and
D = 0 or 1 mod 4), using the precomputed LMFDB fundamental class-group tables.
Falls back to Magma's ClassNumber when the tables do not cover D.}
    require D lt 0 and (D mod 4 in [0, 1]):
        "D must be a negative discriminant (D < 0 and D = 0 or 1 mod 4)";

    order := clGetAssoc(CL_ORDER);
    cached, h := IsDefined(order, D);
    if cached then return h; end if;

    D0 := FundamentalDiscriminant(D);
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
