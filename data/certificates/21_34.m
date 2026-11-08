// Certificates of hyperellipticity for quotients of X_0^21(34).  Read by ReadHyperellipticCertificates.
//
// X_0^21(34)/<w_3, w_14, w_17>, genus 3, shape C: a double cover of the star E = 357d1 (the form is old at 2)
// branched at the rational CM points P_-1428, P_-408, P_-84, P_-4.  Of the three pairings, only
// {P_-1428, P_-84} | {P_-408, P_-4} is certified torsion, by the sharp dual test at level 1428, precision 4000
// (P_-4 has weight 1/4: it is fixed by w_2); it is zero by its class in Phi_3 and Phi_7 and the reduction of
// E(Q)_tors, and the dual graph of Y at 3 is that of the untwisted cover (P_0 = P_-16).  Found 2026-11-08
// (decision 8); #128's model of this curve is not used.
c := AssociativeArray();
c["D"] := 21; c["N"] := 34;
c["W"] := [1, 3, 14, 17, 42, 51, 238, 714];
c["shape"] := "C";
c["E"] := "357d1";
c["divisors"] := [ [<-1428, 1>, <-84, 1>], [<-408, 1>, <-4, 1>] ];
c["dualtests"] := [ <1428, 4000, "sharp", "71f057cdb3a63253525461f61f8d435e"> ];
c["relations"] := [ <1, 2, 1> ];
c["phi"] := [3, 7];
c["tamagawa"] := [3, 7];
c["dualgraph"] := <3, -16>;
c["argument"] := "paper/hyperelliptic-15-38.tex, Sections 1-4 (shape C)";
certificates := [c];
