// Certificates of hyperellipticity for quotients of X_0^15(38).  Read by ReadHyperellipticCertificates.
//
// Y = X_0^15(38)/<w_2, w_15, w_19>, genus 3, shape C: a double cover of the star E = 570e4 branched at the
// rational CM points P_-228, P_-760, P_-12, P_-40.  P_-228 + P_-760 - P_-12 - P_-40 is torsion by the
// sharp dual test at level 1140, precision 4000, on the a_7(E)-eigenspace of T(49), and zero because its
// class in Phi_3 = Z/8 and Phi_5 = Z/4 is zero and the generator of E(Q)_tors = Z/4 and its double leave
// the identity component at 3; the dual graph of Y at 3 is that of the untwisted cover, and the twisted
// cover is not semistable there.  Provenance: PR #151, paper/hyperelliptic-15-38.tex,
// tests/Hyperelliptic15_38.m and tests/_offline/Hyperelliptic15_38.m.
c := AssociativeArray();
c["D"] := 15; c["N"] := 38;
c["W"] := [1, 2, 15, 19, 30, 38, 285, 570];
c["shape"] := "C";
c["E"] := "570e4";
c["divisors"] := [ [<-228, 1>, <-760, 1>], [<-12, 1>, <-40, 1>] ];
c["dualtests"] := [ <1140, 4000, "sharp", "1fbd9ca124791e88a17abc214aebf629", 7> ];
c["relations"] := [ <1, 2, 1> ];
c["phi"] := [3, 5];
c["tamagawa"] := [3, 5];
c["dualgraph"] := <3, -40>;
c["argument"] := "paper/hyperelliptic-15-38.tex, Sections 1-4";
certificates := [c];
