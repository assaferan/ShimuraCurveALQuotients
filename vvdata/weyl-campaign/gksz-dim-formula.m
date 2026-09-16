// Computes dim S_{3/2}(rho_L^*) (equivalently M_{3/2}, cuspidal here -- see
// paper/DRAFT-borcherds-obstruction.md sec.3) via Borcherds' Riemann-Roch formula, GKZ paper
// (Duke 97 (1999), p.9, right after Lemma 4.4):
//   dim HolModForm(rho,k) = d + d*k/12 - alpha(e^{i pi k/2} S) - alpha((e^{i pi k/3} ST)^{-1}) - alpha(T)
// applied to rho = rho_L^* restricted to the Z-compatible ("P=-1", antisymmetric-under-negation)
// eigenspace, at k = 3/2. All quantities are computed via O(n) Gauss sums Sigma_gamma e(c*Q(gamma))
// -- never the full |disc_grp| x |disc_grp| matrix, which is infeasible at this project's scale
// (72,200 to 1,299,272 on the calibration bases below; WeilRepresentationST's matrix build was
// killed after 7+ minutes without finishing the smallest of them).
//
// VALIDATION: all six trace identities below (TS, TS2, TST, TSTST, TS2ST, TS2STST) are checked
// EXACTLY (bit-for-bit, same cyclotomic field) against the real WeilRepresentationST matrices on
// (D,N) = (6,1) [n=72] and (10,1) [n=200] -- see the "for DN in..." loop replaced by the Validate
// calls if you need to re-run that check. Do not trust a further edit to the trace formulas without
// re-running Validate first.
//
// ⚠⚠ RESULT (2026-09-16): this dimension is NOT the deficit.m predictor it was built to be.
// dim M_{3/2}(rho_L^*) is in the THOUSANDS (38_5 -> 1594, 146_1 -> 888) while deficit.m measures
// 1 and 0 there. deficit.m is the RANK of the pairing between this (huge) space and a SMALL fixed
// target set of tracked CM-divisor classes, bounded by dim(target) not dim(S_{3/2}) -- see
// paper/DRAFT-borcherds-obstruction.md sec.5c-5d for the full argument (Waldspurger-type coefficient
// vanishing, not a dimension count -- no closed form is expected to exist for the actual deficit).
// What this file computes is still useful: a real, validated, fast (O(n)) UPPER BOUND on the
// deficit, and the true dimension of the full obstruction space -- not a predictor for the number
// deficit.m returns.

AttachSpec("ShimuraQuotients.spec");

function positive_real_sqrt(K, z, n0, m)
    ee := func<a | z^(Integers()!((a - Floor(a))*n0))>;
    r := K!1; mm := m;
    if IsEven(mm) then
        r *:= ee(1/8) + ee(-1/8);
        mm := mm div 2;
    end if;
    for p in PrimeDivisors(mm) do
        g := &+[ KroneckerSymbol(a, p) * ee(a/p) : a in [1..p-1] ];
        r *:= (p mod 4 eq 1) select g else g * ee(-1/4);
    end for;
    return r;
end function;

// ---- lightweight moments (no matrix build) ----
function RawMoments(Ld)
    D := Ld`D; N := Ld`N;
    dg := Ld`disc_grp; Q := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
    M := IsOdd(D*N) select 4*D*N else 2*D*N;
    sqfree, sq := SquarefreeFactorization(Integers()!Determinant(Q));
    n0 := Lcm([M, 8, 4*sqfree]);
    K<z> := CyclotomicField(n0);
    ee := func<a | z^(Integers()!((a - Floor(a))*n0))>;
    elts := [g : g in dg]; n := #elts;
    vs := [ChangeRing(g@@Ld`to_disc, Rationals()) : g in elts];
    nm := [ (vs[i]*Q, vs[i])/(2*dn^2) : i in [1..n] ];
    invsqrt := 1/(sq * positive_real_sqrt(K, z, n0, sqfree));
    lam := ee(1/4);
    phase := ee(1/8);
    f := 0; torsion_idx := [];
    for i in [1..n] do if 2*elts[i] eq dg!0 then f +:= 1; Append(~torsion_idx, i); end if; end for;
    return n, f, invsqrt, phase, lam, nm, ee, K, n0, torsion_idx;
end function;

function GS(nm, ee, c)
    return &+[ ee(c*x) : x in nm ];
end function;

function RawTraces(Ld)
    n, f, invsqrt, phase, lam, nm, ee, K, n0, torsion_idx := RawMoments(Ld);
    TS      := invsqrt*phase*GS(nm, ee, -2);
    TS2     := lam*f;
    TST     := invsqrt*phase*GS(nm, ee, -1);
    TS2ST   := lam*invsqrt*phase*GS(nm, ee, 3);
    TSTST   := (invsqrt*phase)^2 * GS(nm,ee,1) * GS(nm,ee,-3);
    TS2STST := lam*(invsqrt*phase)^2 * GS(nm,ee,1)^2;
    return n, f, TS, TS2, TST, TS2ST, TSTST, TS2STST, K, n0, invsqrt, phase, lam, nm, ee, torsion_idx;
end function;

// -------------------- validation against real matrices --------------------
procedure Validate(D, N)
    printf "=== validating %o_%o ===\n", D, N;
    Ld := ShimuraCurveLattice(D, N);
    S, T, elts, Km := WeilRepresentationST(Ld);
    n := Nrows(S);
    S2 := S*S; ST := S*T; STST := ST*ST; S2ST := S2*ST; S2STST := S2*STST;
    realTS := Trace(S); realTS2 := Trace(S2); realTST := Trace(ST);
    realTSTST := Trace(STST); realTS2ST := Trace(S2ST); realTS2STST := Trace(S2STST);

    n2, f, TS, TS2, TST, TS2ST, TSTST, TS2STST, K, n0 := RawTraces(Ld);

    printf "n: matrix=%o formula=%o\n", n, n2;
    printf "TS:      real=%o\n         form=%o   equal=%o\n", realTS, TS, realTS eq TS;
    printf "TS2:     real=%o\n         form=%o   equal=%o\n", realTS2, TS2, realTS2 eq TS2;
    printf "TST:     real=%o\n         form=%o   equal=%o\n", realTST, TST, realTST eq TST;
    printf "TSTST:   real=%o\n         form=%o   equal=%o\n", realTSTST, TSTST, realTSTST eq TSTST;
    printf "TS2ST:   real=%o\n         form=%o   equal=%o\n", realTS2ST, TS2ST, realTS2ST eq TS2ST;
    printf "TS2STST: real=%o\n         form=%o   equal=%o\n", realTS2STST, TS2STST, realTS2STST eq TS2STST;
end procedure;

// -------------------- predict dim S_{3/2}(rho_L^*) --------------------
function Predict(D, N : verbose := false)
    Ld := ShimuraCurveLattice(D, N);
    n, f, TS0, TS20, TST0, TS2ST0, TSTST0, TS2STST0, K0, n00, invsqrt, phase, lam0, nm, ee, torsion_idx := RawTraces(Ld);

    // Enlarge to a field containing both the n0-th roots and 12th/8th roots needed for the
    // eigenvalue candidates below.
    n0 := Lcm(n00, 24);
    K<z> := CyclotomicField(n0);
    emb := hom<K0 -> K | z^(n0 div n00)>;
    TS := emb(TS0); TS2 := emb(TS20); TST := emb(TST0);
    TS2ST := emb(TS2ST0); TSTST := emb(TSTST0); TS2STST := emb(TS2STST0);
    lam := emb(lam0);

    conjmap := hom<K -> K | z^(n0-1)>;
    conj := func<x | conjmap(x)>;

    lam_star := conj(lam);
    D_S2 := conj(TS2);
    D_S  := conj(TS);
    D_ST := conj(TST);
    D_STST := conj(TSTST);

    Nminus := (n - D_S2/lam_star)/2;
    M1_S   := (D_S - TS/lam_star)/2;
    M1_ST  := (D_ST - conj(TS2ST)/lam_star)/2;
    M2_ST  := (D_STST - conj(TS2STST)/lam_star)/2;

    // S-block: on the dual's P=-1 eigenspace, S*^2 = -lam_star = ee(1/4), so S* eigenvalues
    // are the square roots a=ee(1/8), b=ee(5/8).
    a := z^(n0 div 8);                // ee(1/8)
    b := z^(5*n0 div 8);              // ee(5/8)
    m_a := (Nminus*b - M1_S)/(b-a);
    m_b := (M1_S - Nminus*a)/(b-a);

    // ST-block: eigenvalues c1=ee(1/12), c5=ee(5/12), c9=ee(9/12); need 12 | n0
    assert n0 mod 12 eq 0;
    c1 := z^(n0 div 12);
    c5 := z^(5*n0 div 12);
    c9 := z^(9*n0 div 12);
    Vm := Matrix(K, 3,3, [1,1,1, c1,c5,c9, c1^2,c5^2,c9^2]);
    rhs := Matrix(K, 3,1, [Nminus, M1_ST, M2_ST]);
    sol := Vm^(-1) * rhs;
    m1 := sol[1,1]; m5 := sol[2,1]; m9 := sol[3,1];

    // twist by e^{i pi k/2} = ee(3/8) at k=3/2: ee(3/8)*a = ee(1/2) (phase 1/2); ee(3/8)*b = ee(1) (phase 0)
    alpha_S  := (K!1/2)*m_a + (K!0)*m_b;
    alpha_ST := (K!2/3)*m1 + (K!1/3)*m5;

    // T-block (dual): eigenvalue phases are frac(-nm[i])
    fullT := 0; torsT := 0;
    for i in [1..n] do
        x := -nm[i]; fr := x - Floor(x);
        fullT +:= fr;
        if i in torsion_idx then torsT +:= fr; end if;
    end for;
    alpha_T := (fullT - torsT)/2;

    dprime := Nminus;
    dim := dprime + dprime*(K!3/2)/12 - alpha_S - alpha_ST - alpha_T;

    if verbose then
        printf "  n=%o f=%o d'=%o\n", n, f, dprime;
        printf "  m_a=%o m_b=%o | m1=%o m5=%o m9=%o\n", m_a,m_b,m1,m5,m9;
        printf "  alpha_S=%o alpha_ST=%o alpha_T=%o\n", alpha_S, alpha_ST, alpha_T;
    end if;
    return dim;
end function;

// Re-validate before trusting any edit above.
Validate(6,1);
Validate(10,1);

bases := [<38,5,1>, <146,1,0>, <194,1,0>, <58,13,2>, <26,31,3>];
for b in bases do
    D:=b[1]; N:=b[2]; known:=b[3];
    printf "\n--- Predict %o_%o (known deficit %o) ---\n", D, N, known;
    t0 := Realtime();
    d := Predict(D, N : verbose := true);
    printf "predicted dim = %o   (%o s)   match=%o (EXPECTED false -- see header)\n",
        d, Realtime(t0), d eq known;
end for;

quit;
