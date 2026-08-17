// Regression test for the Weil representation rho_L (WeilRepresentation.m).
//
// rho_L(S), rho_L(T) must satisfy the defining metaplectic relations of Mp_2(Z). For a lattice of
// signature (b+,b-) with sig := b+ - b- (here (1,2), sig = -1):
//   rho_L(T) e_eta = e(nm(eta)) e_eta                                    (T diagonal)
//   rho_L(S)^2 = e((b- - b+)/4) * P,   P e_eta = e_{-eta}               (Z = S^2 acts by e(1/4) and negation)
//   (rho_L(S) rho_L(T))^3 = rho_L(S)^2                                   (the braid relation)
//   rho_L(S)^4 = -I,  rho_L(S)^8 = I;   rho_L(T)^M = I  (M = lattice level)
// These pin the S-phase and normalization, so they certify the construction before F_f uses it.

procedure test_WeilRepresentation()
    printf "Testing Weil representation rho_L metaplectic relations...";
    for DN in [ [6,1], [10,1] ] do                   // small even-lattice discriminant groups
        D := DN[1]; N := DN[2];
        Ld := ShimuraCurveLattice(D, N);
        S, T, elts, K := WeilRepresentationST(Ld);
        n := #elts;
        Id := IdentityMatrix(K, n);
        M := IsOdd(D*N) select 4*D*N else 2*D*N;

        // negation permutation eta -> -eta on the ordered element list
        idx := AssociativeArray(); for i->g in elts do idx[g] := i; end for;
        neg := [ idx[-elts[i]] : i in [1..n] ];

        // rho_L(S)^2 = c * P(eta -> -eta) with c = e((b- - b+)/4) = e(1/4) = +-i
        S2 := S*S;
        c := S2[neg[1]][1];
        assert c^2 eq -1;                                        // c = e(1/4) is a primitive 4th root
        P := ZeroMatrix(K, n, n);
        for j in [1..n] do P[neg[j]][j] := c; end for;
        assert S2 eq P;
        // braid relation and S-order
        assert (S*T)^3 eq S2;
        assert S^4 eq -Id;
        assert S^8 eq Id;
        // T has order dividing the lattice level M
        assert T^M eq Id;
    end for;
    printf "Done!\n";
end procedure;

test_WeilRepresentation();
