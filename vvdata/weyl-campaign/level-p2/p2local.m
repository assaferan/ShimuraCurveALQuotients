// Local structure of the trace-zero lattice of an Eichler order at a prime p dividing the level,
// with p || N as a control and p^2 || N as the case Yang's Lemma 18 does not cover.
AttachSpec("ShimuraQuotients.spec");
for b in [[6,5,5],[6,25,5],[6,125,5],[15,2,2],[15,4,2],[15,8,2],[14,3,3],[14,9,3],[10,9,3],[21,4,2]] do
    D := b[1]; N := b[2]; p := b[3];
    B := QuaternionAlgebra(D); O := QuaternionOrder(B, N);
    bas := Basis(O);
    Mz := Matrix(Integers(), 4, 1, [Integers() | Trace(x) : x in bas]);
    K := KernelMatrix(Mz);
    L0 := [&+[K[i][j]*bas[j] : j in [1..4]] : i in [1..Nrows(K)]];
    gram := Matrix(Integers(), 3, 3, [Integers() | Trace(L0[i]*Conjugate(L0[j])) : i, j in [1..3]]);
    ed := ElementaryDivisors(gram);
    printf "D=%-3o N=%-4o p=%o : reduced disc %o, det(Gram) = %o, p-valuations of the elementary divisors %o\n",
           D, N, p, Discriminant(O), Factorization(Determinant(gram)), [Valuation(e, p) : e in ed];
end for;
