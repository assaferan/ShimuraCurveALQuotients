// Certificates for quotients of X_0^141(2).  Read by ReadHyperellipticCertificates.
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
certificates := [c];
