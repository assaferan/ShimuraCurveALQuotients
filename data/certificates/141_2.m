// Certificates for quotients of X_0^141(2): one not hyperelliptic (shape N), one hyperelliptic (shape C).  Read by ReadHyperellipticCertificates.
//
// X_0^141(2)/{1, 6, 94, 141}, genus 3, shape N (NOT hyperelliptic): a double cover of the genus-one star
// branched at the rational CM points P_-4, P_-12 and the place Z_-1128 of degree two, so the only
// Galois-stable pairing is {P_-4, P_-12} | Z_-1128; its difference has class 2 in Phi_3 = Z/7, so it is not
// zero and the hyperelliptic involution, which would give such a pairing, does not exist.  This is the 34th
// of paper/hyperelliptic-candidates.tex; #128 had decided it by the enumeration of the Heegner-height
// route, which this certificate does not use.
c := AssociativeArray();
c["D"] := 141; c["N"] := 2;
c["W"] := [1, 6, 94, 141];
c["shape"] := "N";
c["phi"] := [3];
c["argument"] := "paper/hyperelliptic-15-38.tex, Section 1 (the criterion); paper/component-groups.tex";
// X_0^141(2)/{1, 3, 94, 282}, genus 3, shape C (hyperelliptic): a double cover of the star E = 141a1 (the
// form is old at 2) branched at P_-24, P_-4 (rational) and Z_-564 (degree two); P_-24 + P_-4 - Z_-564 is
// torsion by the sharp dual test at level 564, precision 4000 (P_-4 has weight 1/4: it is fixed by w_2), zero
// by its class in Phi_3 and Phi_47 and the reduction of E(Q)_tors, and the dual graph of Y at 3 is that of
// the untwisted cover (P_0 = P_-7).  Found 2026-11-08 (decision 8); #128's model is not used.
h := AssociativeArray();
h["D"] := 141; h["N"] := 2;
h["W"] := [1, 3, 94, 282];
h["shape"] := "C";
h["E"] := "141a1";
h["divisors"] := [ [<-24, 1>, <-4, 1>], [<-564, 1>] ];
h["dualtests"] := [ <564, 4000, "sharp", "dae4b81d4be5314adb9ce54681620d77"> ];
h["relations"] := [ <1, 2, 1> ];
h["phi"] := [3, 47];
h["tamagawa"] := [3, 47];
h["dualgraph"] := <3, -7>;
h["argument"] := "paper/hyperelliptic-15-38.tex, Sections 1-4 (shape C)";
certificates := [c, h];
