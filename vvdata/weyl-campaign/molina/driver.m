// Driver: Molina's Steps 1-4 end to end for one curve, with verification.
//   magma -b D:=39 m:=1 g:=3 logfile:=out.txt driver.m < /dev/null
// Optional: maxj2 (bound on the 2-adic part of each J_i, default 2); maxidx (skip an
// IndexFormEquation above this index); onlyfirst (use only the smallest |a_R|); dinf; verbose.

AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
load "vvdata/weyl-campaign/molina/molina.m";

D := StringToInteger(D);
m := StringToInteger(m);
g := StringToInteger(g);
maxj2 := assigned maxj2 select StringToInteger(maxj2) else 2;
verbose := assigned verbose select StringToInteger(verbose) else 1;
forced := assigned dinf select StringToInteger(dinf) else 0;
maxidx := assigned maxidx select StringToInteger(maxidx) else 200000;
onlyfirst := assigned onlyfirst select StringToInteger(onlyfirst) else 1;

// Magma buffers stdout when it is a file, so progress goes to its own handle, flushed per line.
LOG := Open(logfile, "w");
procedure Say(s)
    fprintf LOG, "%o\n", s;
    Flush(LOG);
end procedure;

W := (m eq 1) select {Integers() | 1} else {Integers() | 1, m};
Say(Sprintf("=== X_0^%o(1)%o, genus %o", D, (m eq 1) select "" else Sprintf("/<w_%o>", m), g));
assert GenusShimuraCurveQuotient(D, 1, W) eq g;

// ---------------- Step 1: the Weierstrass decomposition
wd := MolWeierstrassData(D, m, g);
if #wd ne 1 then
    // Molina's Section 10.2 assumes the hyperelliptic involution of Y is the image of some w_n;
    // when no divisor n of D makes the fixed points add up to 2g+2, the method does not apply.
    Say(Sprintf("NOT APPLICABLE: %o divisors of D give a set of Weierstrass points of size 2g+2",
        #wd));
    Say("DRIVER DONE");
    quit;
end if;
n := wd[1][1]; np := wd[1][2]; rows := wd[1][3];
// The two points at infinity are pi(P_inf) and pi(w_n(P_inf)), and they have to be exchanged by
// complex conjugation, which acts on Heegner points as w_D.  That needs w_D to agree with one of
// 1, w_m, w_n, w_n' on CM(R_inf), i.e. D in {m, n, n'}; otherwise the pair is not Galois stable
// and Molina's choice of points at infinity (Sections 3 and 10.2) is unavailable.
if not (D in {m, n, np}) then
    Say(Sprintf("NOT APPLICABLE: the hyperelliptic involution is w_%o and D is not one of %o,"
        * " so the two points at infinity are not exchanged by complex conjugation", n, {m, n, np}));
    Say("DRIVER DONE");
    quit;
end if;
Say(Sprintf("Step 1: hyperelliptic involution w_%o (paired with w_%o); orbits %o", n, np, rows));

X := CreateShimuraQuot(D, 1, W);
facs := [* *];
rows := Sort(rows, func<a, b | b[2] - a[2]>);     // the largest factor anchors the assembly
for r in rows do
    d := r[1]; deg := r[2]; mult := r[3];
    if deg eq 1 then
        F := Rationals(); OF := Integers(); dF := 1;
    else
        // the ring class field machinery hands back a presentation with huge coefficients,
        // which makes IndexFormEquation unusable
        F := OptimizedRepresentation(AbsoluteField(FieldOfDefinitionOfCMPoint(X, d)));
        OF := MaximalOrder(F);
        dF := Discriminant(OF);
    end if;
    for i in [1 .. mult] do Append(~facs, <d, deg, F, OF, dF>); end for;
end for;
Say(Sprintf("  R(x) = %o factors of degrees %o over fields of discriminant %o",
    #facs, [f[2] : f in facs], [f[5] : f in facs]));
assert &+[f[2] : f in facs] eq 2*g+2;

vdisc, Dtarget, th := MolDiscR(D, m, g);
Say(Sprintf("Step 1: |disc R| = %o", Factorization(Dtarget)));

// ---------------- Step 2: the leading coefficient
Sprimes := Sort([2] cat PrimeDivisors(D));
aRs := [];
for dd in MolInfinityDiscs(D) do
    if forced ne 0 and dd ne forced then continue; end if;
    aR, ok, sq := MolLeadingCoefficient(D, m, n, np, dd);
    Say(Sprintf("Step 2: K_inf of discriminant %o -> a_R = %o (square class %o, consistent: %o)",
        dd, aR, sq, ok));
    if ok then Append(~aRs, aR); end if;
end for;
aRs := Sort(SetToSequence(Seqset(aRs)), func<a,b | Abs(a) - Abs(b)>);
assert #aRs gt 0;
if onlyfirst eq 1 then aRs := [aRs[1]]; end if;
Say(Sprintf("  using a_R in %o", aRs));

// ---------------- enumerations
function AbsSplits(N, k)
    if k eq 1 then return [[N]]; end if;
    out := [];
    for dv in Divisors(N) do
        for rest in $$(N div dv, k-1) do Append(~out, [dv] cat rest); end for;
    end for;
    return out;
end function;

function SignPatterns(k)
    if k eq 1 then return [[1], [-1]]; end if;
    out := [];
    for rest in $$(k-1) do Append(~out, [1] cat rest); Append(~out, [-1] cat rest); end for;
    return out;
end function;

function SignedSplits(aR, k)
    out := [];
    for v in AbsSplits(Abs(aR), k) do
        for s in SignPatterns(k) do
            sv := [s[i]*v[i] : i in [1..k]];
            if &*sv eq aR then Append(~out, sv); end if;
        end for;
    end for;
    return Sort(SetToSequence(Seqset(out)));
end function;

function Tuples(pools)
    if #pools eq 0 then return [[]]; end if;
    out := [];
    for rest in $$(pools[2..#pools]) do
        for v in pools[1] do Append(~out, [v] cat rest); end for;
    end for;
    return out;
end function;

function JTuples(facs, vdisc, Dtarget, Sprimes, maxj2, th)
    pools := [];
    for f in facs do
        if f[2] eq 1 then Append(~pools, [1]); continue; end if;
        dM := Abs(f[5]);
        exps := [];
        for p in Sprimes do
            hi := (Min(vdisc[p], MolDiscBound(th, p, f[2], vdisc[p])) - Valuation(dM, p)) div 2;
            if p eq 2 then hi := Min(hi, maxj2); end if;
            if hi lt 0 then return []; end if;
            Append(~exps, hi);
        end for;
        pool := [&*[Sprimes[i]^t[i] : i in [1..#Sprimes]] : t in Tuples([[0..e] : e in exps])];
        Append(~pools, Sort(pool));
    end for;
    out := [];
    for J in Tuples(pools) do
        prod := &*[Integers() | J[i]^2 * Abs(facs[i][5]) : i in [1..#facs]];
        if Dtarget mod prod ne 0 then continue; end if;
        if not IsSquare(Dtarget div prod) then continue; end if;
        Append(~out, J);
    end for;
    return out;
end function;

// ---------------- Step 4
cache := AssociativeArray();
procedure Cands(f, a, J, cap, ~cache, ~out)
    key := <f[5], f[2], a, J>;
    if IsDefined(cache, key) then out := cache[key]; return; end if;
    I := MolIndex(a, f[2], J);
    if f[2] ge 2 and I gt cap then
        Say(Sprintf("    SKIPPED index %o (field disc %o) -- above the cap", I, f[5]));
        cache[key] := [];
        out := [];
        return;
    end if;
    t := Cputime();
    c := MolFactorCandidates(f[3], f[4], a, f[2], J);
    if f[2] ge 2 then
        Say(Sprintf("    IndexFormEquation(field disc %o, index %o) -> %o candidates, %o s",
            f[5], I, #c, Round(Cputime(t))));
    end if;
    cache[key] := c;
    out := c;
end procedure;

Xq := CreateShimuraQuot(D, 1, W);
Xq`g := g;
found := {Qx | };
bnd := Isqrt(Dtarget) + 1;

function Verify(R)
    if Degree(R) ne 2*g+2 then return false; end if;
    if Abs(Integers() ! Discriminant(R)) ne Dtarget then return false; end if;
    C := HyperellipticCurve(Qx ! R);
    if Genus(C) ne g then return false; end if;
    ntest := 0;
    for p in [q : q in PrimesUpTo(30) | q ne 2 and D mod q ne 0] do
        fp := ChangeRing(Qx ! R, GF(p));
        if Discriminant(fp) eq 0 or Degree(fp) lt 2*g+1 then continue; end if;
        if #Points(ChangeRing(C, GF(p))) ne ComputePointsViaTrace(Xq, p, 1) then return false; end if;
        if #Points(ChangeRing(C, GF(p^2))) ne ComputePointsViaTrace(Xq, p, 2) then return false; end if;
        ntest +:= 1;
    end for;
    return ntest ge 2;
end function;

// disc(R) = prod disc(p_i) . prod Res(p_i,p_j)^2, and both disc R and the disc(p_i) = J_i^2 disc(M_i)
// are known, so the product of the pairwise resultants is known exactly.  Each one therefore
// divides it -- a far tighter test than "supported on 2 and the primes of D".
function GoodRes(r, Rtot)
    if r eq 0 then return false; end if;
    return Rtot mod Abs(r) eq 0;
end function;

nprod := 0;
ncombo := 0;
Js := JTuples(facs, vdisc, Dtarget, Sprimes, maxj2, th);
Say(Sprintf("Step 3: %o admissible discriminant splittings (J-tuples)", #Js));

for aR in aRs do
  spl := SignedSplits(aR, #facs);
  Say(Sprintf("a_R = %o: %o splittings of the leading coefficient", aR, #spl));
  for av in spl do
    for J in Js do
        lists := [];
        ok := true;
        for i in [1 .. #facs] do
            c := [];
            Cands(facs[i], av[i], J[i], maxidx, ~cache, ~c);
            if #c eq 0 then ok := false; break; end if;
            Append(~lists, c);
        end for;
        if not ok then continue; end if;
        dprod := &*[Integers() | J[i]^2 * Abs(facs[i][5]) : i in [1 .. #facs]];
        Rtot := Isqrt(Dtarget div dprod);
        ncombo +:= 1;
        if verbose ge 2 then
            Say(Sprintf("  a = %o, J = %o, candidates %o, prod|Res| = %o",
                av, J, [#l : l in lists], Rtot));
        end if;
        for c1 in lists[1] do
            stack := [* <[c1], 2> *];
            while #stack gt 0 do
                st := stack[#stack]; Prune(~stack);
                cur := st[1]; idx := st[2];
                if idx gt #facs then
                    nprod +:= 1;
                    R := &*cur;
                    if Verify(R) then
                        Include(~found, Qx ! R);
                        Say(Sprintf("  *** y^2 = %o", R));
                    end if;
                    continue;
                end if;
                for c in lists[idx] do
                    for r in MolTranslationRange(cur[1], c, Rtot) do
                        cc := Evaluate(c, x + r);
                        if &and[GoodRes(Integers() ! Resultant(u, cc), Rtot) : u in cur] then
                            Append(~stack, <cur cat [cc], idx+1>);
                        end if;
                    end for;
                end for;
            end while;
        end for;
    end for;
  end for;
end for;

Say(Sprintf("RESULT D=%o m=%o: %o equation(s); %o (a,J) combinations, %o products assembled",
    D, m, #found, ncombo, nprod));
for R in found do Say(Sprintf("  y^2 = %o", R)); end for;
Say("DRIVER DONE");
quit;
