// tests/HyperellipticCRV.m
//
// Which of the curves stored as a list of equations (the "CRV" shape: fibre products of double
// covers of the star line, and the older pairs over a conic) are hyperelliptic, and whether any of
// them has a model y^2 = f(x) over Q that the file should hold instead.
//
// The test is Max Noether's count of quadrics through the canonical image (IsHyperellipticFunctionField):
// (g-1)(g-2)/2 linear relations among the pairwise products of the holomorphic differentials for a
// hyperelliptic curve, (g-2)(g-3)/2 otherwise.  For every hyperelliptic curve found, the test then
// checks two things.  Its hyperelliptic involution is an Atkin-Lehner one: the product of the two
// double covers, or one of them, is, up to squares, of degree at most 2, so the curve is the double
// cover of the genus-0 quotient y^2 = c(t) by one of its three involutions (a hyperelliptic
// involution outside that group would be new, and the test fails so that it is looked at).  And
// that conic has no rational point: if it had one, the curve would have a model y^2 = f(x) over Q
// (HyperellipticModelOverConic), and the file should store that instead -- as the two covers of
// X_0(10,53) by w_530 and w_53 now do.  A hyperelliptic curve of even genus always has such a model.
//
// Both checks assume exactly two factors.  With k factors the double covers y^2 = f_S(t) are the
// quotients by the index-2 subgroups of (Z/2)^k, and the quotient by an involution w, the one that
// decides whether the hyperelliptic involution is Atkin-Lehner and whose rational points decide the
// model, is among them only when k = 2; for k >= 3 it has degree 2^(k-1) over the line.  So the test
// fails if a hyperelliptic entry has three or more factors, rather than conclude anything from it.
//
// Controls: y^2 = h(t), w^2 = t, z^2 = t - 1 with deg h = 4 (genus 7, the hyperelliptic curve
// y^2 = h(w^2) together with z) must be detected, and the three plane quartics X_0(6,23)/w_23,
// X_0(34,3)/w_2, X_0(46,3)/w_3 (data/models/PROVENANCE.md, 2026-10-03) must not be.
//
// Sweep of 2026-10-14 (main at a2e0747 with the two 10_53 covers converted): 297 CRV entries,
// 58 hyperelliptic, every one by an Atkin-Lehner involution over a conic without a rational point.

printf "Hyperellipticity of every stored CRV curve...";
_ := ClassNumberLU(-4);
P<x> := PolynomialRing(Rationals());

// the controls
hyp_ctrl := IsHyperellipticFunctionField(FibreProductFunctionField([x^4 + 3*x^3 - x + 5, x, x - 1]));
assert hyp_ctrl;

crv_files := Split(Pipe("ls data/models/*.m | xargs grep -l '\"CRV\"'", ""), "\n");
n_crv := 0; n_hyp := 0; bad := []; quartics_seen := 0;
quartics := {<6, 23, {Integers()|1, 23}>, <34, 3, {Integers()|1, 2}>, <46, 3, {Integers()|1, 3}>};
for f in crv_files do
    base := Split(Split(f, "/")[3], ".")[1]; base := base[8..#base];
    D := StringToInteger(Split(base, "_")[1]); N := StringToInteger(Split(base, "_")[2]);
    models := eval (Read(f) cat "\nreturn models;");
    for k in Keys(models) do
        for e in models[k] do
            if Type(e[2]) ne MonStgElt then continue; end if;
            n_crv +:= 1;
            fs := CRVDoubleCovers(e[3]);
            K := FibreProductFunctionField(fs);
            g := Genus(K);
            assert g eq e[1];
            hyp := IsHyperellipticFunctionField(K);
            if <D, N, Set(k)> in quartics then
                quartics_seen +:= 1;
                assert not hyp;
            end if;
            if not hyp then continue; end if;
            n_hyp +:= 1;
            tag := Sprintf("%o %o (genus %o)", base, Sort(SetToSequence(Set(k))), g);
            if #fs ne 2 then
                Append(~bad, tag cat ": hyperelliptic with " cat Sprint(#fs) cat " factors; the conic test decides the involution only for two");
                continue;
            end if;
            ok, c := ConicFactorOfFibreProduct(fs);
            if not ok then
                Append(~bad, tag cat ": hyperelliptic, but by an involution outside the Atkin-Lehner group");
            elif IsEven(g) or ConicHasRationalPoint(c) then
                Append(~bad, tag cat ": a double cover of a conic with a rational point, so it has a model y^2 = f(x) over Q");
            end if;
        end for;
    end for;
end for;

error if n_crv eq 0, "HyperellipticCRV: no CRV entries found, nothing was checked";
error if quartics_seen ne 3, "HyperellipticCRV: the three plane quartics were not all found";
error if not IsEmpty(bad),
    Sprintf("HyperellipticCRV: %o stored curve(s) should be stored differently: %o", #bad, bad);
printf " ok (%o CRV entries, %o hyperelliptic, each with two factors and a double cover of a conic without a rational point)\n", n_crv, n_hyp;
