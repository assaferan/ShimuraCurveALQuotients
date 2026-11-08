// Certificates of hyperellipticity for quotients of X_0^6(341).  Read by ReadHyperellipticCertificates.
//
// Y = X_0^6(341)/<w_3, w_11, w_62>, genus 4, shape B2: a double cover of the star (genus 2) branched at
// P_-24 + P_-88, related to P_-43 + P_-123 and P_-84 + P_-228 by three Borcherds forms; the star's
// Jacobian J has rational torsion of order at most 2, so the branch pair is K or K + t.  The Hecke algebra
// acts on J through Z[sqrt 2] (T_23 = -4 +- sqrt 2), and the two principally polarised surfaces with that
// real multiplication in the isogeny class are the Jacobians of the curves below, C (the star, fitted from
// its differentials) and C_A = Jac(C)/<t>; on neither is K + t twice a rational class.  The real
// multiplication is CHIMP's certified endomorphism (tangent matrices below, logs in tools/proof-6_341/
// step2-endomorphisms*.log), the conductor exponent at 2 CHIMP's proved one (step0-conductor-at-2*.log).
c := AssociativeArray();
c["D"] := 6; c["N"] := 341;
c["W"] := [1, 3, 11, 33, 62, 186, 682, 2046];
c["shape"] := "B2";
c["divisors"] := [ [<-24, 1>, <-88, 1>], [<-43, 1>, <-123, 1>], [<-84, 1>, <-228, 1>] ];
c["forms"] := "data/borcherds/forms_6_341.m";
c["torsion"] := [60, 6];
c["heckeprime"] := 23;
// y^2 + h y = f, coefficients from the constant term
c["stars"] := [ <[9, -33, 9, 29, -12, 6], [0, 1, 1]>, <[144, -120, 265, -82, 133, -12, 36], []> ];
c["rm"] := [ [0, 1, 2, 0], [-1, 1, 1, 1] ];
c["conductor2"] := [2, 2];
c["reducibility"] := [4093, 12277, 14323, 16369];
c["argument"] := "paper/hyperelliptic-6-341.tex; paper/classification-certificates.tex";
certificates := [c];
