// Certificates of hyperellipticity for quotients of X_0^158(3).  Read by ReadHyperellipticCertificates.
//
// X_0^158(3)/<w_2, w_237>, genus 3, shape O: a double cover of the star, whose Jacobian is 158a1 (the form is
// old at 3), branched at the four CM points of discriminant -632, one place over a D4 field with quadratic
// subfield Q(sqrt -79).  3 splits in Q(sqrt -158) and chi(p_3) = (-79/3) = -1 for the genus character of
// Q(sqrt -79), so the chi-component of the branch divisor vanishes on the star and the stable pairing has
// equal sums; 158a1 and its twist by -79 have no torsion.  Not one of the 34: the special fibre at 3 cannot rule
// it out (Z = X_0^158(1)/<w_2> is not hyperelliptic over Q, where the same component is nonzero).  Found
// 2026-11-08.
c := AssociativeArray();
c["D"] := 158; c["N"] := 3;
c["W"] := [1, 2, 237, 474];
c["shape"] := "O";
c["E"] := "158a1";
c["p"] := 3;
c["ratpoint"] := -8;
c["argument"] := "paper/classification-certificates.tex, Section 6 (shape O)";
certificates := [c];
