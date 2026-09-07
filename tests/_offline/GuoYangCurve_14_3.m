// tests/_offline/GuoYangCurve_14_3.m
//
// Does our X_0^14(3) model reproduce Guo-Yang's published curve?
//
//     Guo-Yang:  z^2 = -9x^2 - 2  ,  y^2 = -7x^4 + 22x^2 + 1        [genus 3]
//
// ⚠ WHY THIS IS OFFLINE AND NOT IN tests/GuoYangEquations.m. The comparison is a single
// IsIsomorphic on a genus-3 complete intersection in weighted projective space, and it takes
// **6817 s (1 h 54 m)**. GuoYangEquations.m runs in ~97 s for nine bases; adding this would make
// it unusable in CI. Cost here is driven by the PRESENTATION, not the genus -- 39_2's genus-7
// W={1} is plain hyperelliptic and settles in 0.06 s, while 21_2's genus-3 CRV pair takes ~100 s
// and this one takes two hours.
//
// PROVENANCE OF THE MODEL. data/models/models_14_3.m was produced with CMNONCOPRIME=1; under
// default settings 14_3's covers are under-determined and W={1} comes out EMPTY. The flag has no
// theoretical guarantee (the p | gcd(d,N) local factor has no live implementation), so this
// comparison IS the justification for committing that file. If the CM filter changes, run this.
//
// ⚠⚠ THE TRAP THIS TEST EXISTS TO OUTLIVE. Before the exact test finished I tried to settle the
// question cheaply by comparing quotients, and concluded the curves were NOT isomorphic. That was
// WRONG. For y^2=f (deg 4), x^2=g (deg 2) the three double-cover quotients have genus 0, 1, 2 and
// it is tempting to infer that an isomorphism must match ours to theirs. That needs the Klein
// 4-group to be UNIQUE in Aut(C), and here it is not: D*N = 42, so W={1} has SEVEN involution
// quotients (genera 1,2,2,2,0,1,1) and several valid V4s. Measured: our W={1} is presented over
// V4 = {1,2,7,14}, Guo-Yang's over {1,3,14,42}. Comparing presentation-to-presentation gives
// false on BOTH quotients, and the degree-4 branch forms even have different Galois groups
// (C2^2 vs D4) -- all of which refutes only a BASE-COMPATIBLE isomorphism for that pairing.
// The right question is whether each of THEIR quotients matches ANY of ours:
//       GY genus-1 quotient ~ our [1,42]        GY genus-2 quotient ~ our [1,3]
// which is consistent -- and the exact test below then returned true.
//
// Nor does a passing point-count screen settle it: tests/IsoScreen.m gives 13/13 consistent here,
// but non-isomorphic curves with isogenous Jacobians would do the same. Only this is proof.
//
// Run:  magma -b filename:=tests/_offline/GuoYangCurve_14_3.m run_tests.m    (~2 h)

printf "X_0^14(3): comparing our W={1} model to Guo-Yang's published curve (SLOW, ~2h)...\n";

gy14_models := eval (Read("data/models/models_14_3.m") cat "\nreturn models;");
gy14_ok, gy14_entry := IsDefined(gy14_models, [Integers()|1]);
error if not gy14_ok, "models_14_3.m has no W={1} key";
error if #gy14_entry eq 0,
    "models_14_3.m has an EMPTY W={1} entry -- that is what the default CM filter produces; "
    * "this file needs the CMNONCOPRIME=1 model (see the model file header)";
gy14_g, gy14_tag, gy14_strs := Explode(gy14_entry[1]);
error if Type(gy14_tag) ne MonStgElt, "expected a stored CRV pair for W={1}";

// y has weight 2: its equation is degree 4.
Pours14<x,y,z,s> := WeightedProjectiveSpace(Rationals(), [1,2,1,1]);
gy14_Cours := Curve(Pours14, [ eval ("return " cat str cat ";") : str in gy14_strs ]);

Q14<a,b,c,d> := WeightedProjectiveSpace(Rationals(), [1,2,1,1]);
// a = x, b = y, c = z, d = homogenising variable
gy14_Cgy := Curve(Q14, [ b^2 - (-7*a^4 + 22*a^2*d^2 + d^4),
                         c^2 - (-9*a^2 - 2*d^2) ]);

error if Genus(gy14_Cours) ne Genus(gy14_Cgy),
    Sprintf("X_0^14(3): our genus %o vs Guo-Yang's %o -- comparing the wrong object",
            Genus(gy14_Cours), Genus(gy14_Cgy));

gy14_t0 := Cputime();
error if not IsIsomorphic(gy14_Cours, gy14_Cgy),
    "X_0^14(3): our W={1} model is NOT isomorphic to Guo-Yang's published curve";
printf "  ok -- genus %o, IsIsomorphic in %o s\n", Genus(gy14_Cours), Cputime(gy14_t0);
