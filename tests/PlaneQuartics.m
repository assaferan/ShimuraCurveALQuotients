// tests/PlaneQuartics.m
//
// THE THREE QUOTIENTS RECORDED AS NOT HYPERELLIPTIC: X_0^6(23)/w_23, X_0^34(3)/w_2, X_0^46(3)/w_3.
//
// Each is a genus-3 quotient all of whose three index-2 Atkin-Lehner quotients have genus 1, so it
// is the fibre product over the star line of any two of them, and the model files store those
// elliptic quotients as y^2 + h y = f.  The verdict "not hyperelliptic" recorded for these three
// in data/models/PROVENANCE.md and HANDOFF.md (2026-10-03) is why they are stored as fibre
// products and have no equation y^2 = f(x); nothing in the suite checked it.
//
// For each of the three this test rebuilds the curve from the stored elliptic quotients
// (tools/planequartic.m is the same computation run by hand) and checks four things: the genus is
// the one the Shimura-curve genus formula gives; the point counts over F_p and F_{p^2} at the
// first two good primes agree with the Eichler-Selberg trace formula on the W-fixed part of the
// D-new space (ComputePointsViaTrace, which is evaluated from D, N and W and knows nothing about
// the model); the curve is not hyperelliptic; and its canonical image is a nonsingular plane
// quartic, which is the geometric content of the verdict -- a hyperelliptic genus-3 curve has a
// conic for its canonical image.  A prime at which the tower cannot be built or the genus drops is
// bad reduction of the model and is replaced by the next, searching the first 40.
//
// The hyperellipticity check is not vacuous: run on X_0^6(23)/w_2, the genus-3 hyperelliptic
// curve stored in the same file, the same two calls return "hyperelliptic" and a canonical image
// of degree 2 in P^2, a conic (control, 2026-10-13).

printf "Checking the three quotients recorded as not hyperelliptic...\n";

pq_targets := [<6, 23, 23>, <34, 3, 2>, <46, 3, 3>];
pq_curves := GetHyperellipticCandidates();
pq_nprime := 0;

for pq_T in pq_targets do
    pq_D := pq_T[1]; pq_N := pq_T[2]; pq_w := pq_T[3];
    pq_W := {Integers()| 1, pq_w};
    pq_tag := Sprintf("X_0^%o(%o)/w_%o", pq_D, pq_N, pq_w);
    P<x> := PolynomialRing(Rationals());
    pq_models := eval (Read(Sprintf("data/models/models_%o_%o.m", pq_D, pq_N)) cat "\nreturn models;");
    // the stored double covers of the star curve above W, all of genus 1 here
    pq_poly := AssociativeArray(); pq_avail := {};
    for pq_U in AtkinLehnerDoubleCoversOver(pq_W, pq_D, pq_N) do
        pq_k := Sort(SetToSequence(pq_U));
        if not IsDefined(pq_models, pq_k) or #pq_models[pq_k] eq 0
           or Type(pq_models[pq_k][1][2]) eq MonStgElt then continue; end if;
        pq_e := pq_models[pq_k][1];
        error if pq_e[1] ne 1,
            Sprintf("PlaneQuartics: %o has a quotient of genus %o above it, not 1", pq_tag, pq_e[1]);
        pq_poly[pq_U] := pq_e[2] + (#pq_e ge 3 select pq_e[3]^2/4 else 0);
        Include(~pq_avail, pq_U);
    end for;
    pq_ok, pq_gens := FibreProductGenerators(pq_W, pq_avail, pq_D, pq_N);
    error if not pq_ok,
        Sprintf("PlaneQuartics: the %o stored elliptic quotients above %o do not generate it",
                #pq_avail, pq_tag);
    pq_fs := [pq_poly[pq_U] : pq_U in pq_gens];
    error if #pq_fs ne 2,
        Sprintf("PlaneQuartics: %o needs %o factors, not two", pq_tag, #pq_fs);

    pq_A := AffineSpace(Rationals(), 3);
    pq_C := ProjectiveClosure(Curve(pq_A, [pq_A.(1+i)^2 - Evaluate(pq_fs[i], pq_A.1) : i in [1..2]]));
    pq_g := Genus(pq_C);
    pq_X := rep{Y : Y in pq_curves | Y`D eq pq_D and Y`N eq pq_N and Y`W eq pq_W};
    error if pq_g ne pq_X`g or pq_g ne 3,
        Sprintf("PlaneQuartics: %o has genus %o, the genus formula gives %o", pq_tag, pq_g, pq_X`g);

    pq_done := 0; pq_p := 2; pq_tried := 0;
    while pq_done lt 2 and pq_tried lt 40 do
        pq_p := NextPrime(pq_p); pq_tried +:= 1;
        if (pq_D * pq_N) mod pq_p eq 0 then continue; end if;
        pq_Kp := RationalFunctionField(GF(pq_p)); pq_L := pq_Kp; pq_built := true;
        try
            for pq_f in pq_fs do
                RY<Y> := PolynomialRing(pq_L);
                pq_L := FunctionField(Y^2 - pq_L!Evaluate(PolynomialRing(GF(pq_p))!pq_f, pq_Kp.1));
            end for;
        catch err pq_built := false;
        end try;
        if not pq_built or Genus(pq_L) ne pq_g then continue; end if;
        pq_cnt := [&+[e * #Places(pq_L, e) : e in Divisors(d)] : d in [1..2]];
        pq_exp := [ComputePointsViaTrace(pq_X, pq_p, d) : d in [1..2]];
        error if pq_cnt ne pq_exp,
            Sprintf("PlaneQuartics: %o over F_%o, F_%o^2 has %o points, the trace formula gives %o",
                    pq_tag, pq_p, pq_p, pq_cnt, pq_exp);
        pq_done +:= 1; pq_nprime +:= 1;
    end while;
    error if pq_done lt 2,
        Sprintf("PlaneQuartics: %o had only %o good primes among the first 40", pq_tag, pq_done);

    error if IsHyperelliptic(pq_C),
        Sprintf("PlaneQuartics: %o is hyperelliptic, against the verdict recorded in PROVENANCE.md",
                pq_tag);
    pq_img := CanonicalImage(pq_C, CanonicalMap(pq_C));
    error if Dimension(Ambient(pq_img)) ne 2 or Degree(pq_img) ne 4 or not IsNonsingular(pq_img),
        Sprintf("PlaneQuartics: the canonical image of %o is not a nonsingular plane quartic: "
                * "it lies in P^%o and has degree %o", pq_tag,
                Dimension(Ambient(pq_img)), Degree(pq_img));
    printf "  %o: genus 3, not hyperelliptic, canonical image a nonsingular plane quartic\n", pq_tag;
end for;

error if pq_nprime ne 2 * #pq_targets,
    Sprintf("PlaneQuartics: NO EVIDENCE -- %o point-count comparisons, expected %o",
            pq_nprime, 2 * #pq_targets);

printf "  ok (%o curves, %o primes each over F_p and F_p^2)\n", #pq_targets, 2;
