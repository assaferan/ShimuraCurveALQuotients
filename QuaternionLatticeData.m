declare type QuaternionLatticeData;

declare attributes QuaternionLatticeData: D, N, Q, O, L, Ldual, Qinv, basis_L, disc_grp, to_disc;

function get_lattice_data(D, N)
    B := QuaternionAlgebra(D);
    O_max := MaximalOrder(B);
    O := Order(O_max,N);
    basis_O := Basis(O);
    L_space := Kernel(Transpose(Matrix(Integers(),[[Trace(x) : x in basis_O]])));
    basis_L := [&+[b[i]*basis_O[i] : i in [1..4]] : b in Basis(L_space)];
    BM_L := Matrix([Eltseq(b) : b in basis_L]);
    Q := Matrix([[Norm(x+y)-Norm(x)-Norm(y) : y in basis_L] : x in basis_L]);
    BM_Ldual := Q^(-1)*BM_L;
    // L := LatticeWithGram(Q : CheckPositive := false);
    // return L;
    denom := Denominator(BM_Ldual);
    // We are modifying it to be always with respect to the basis of L.
    // Ldual := RSpaceWithBasis(ChangeRing(denom*BM_Ldual,Integers()));
    Ldual := RSpaceWithBasis(ChangeRing(denom*Q^(-1), Integers()));
    // L := RSpaceWithBasis(ChangeRing(denom*BM_L,Integers()));
    L := RSpaceWithBasis(ScalarMatrix(3,denom));
    disc_grp, to_disc := Ldual / L;
    return L, Ldual, disc_grp, to_disc, Q^(-1), Q, O, basis_L;
end function;

intrinsic ShimuraCurveLattice(D::RngIntElt, N::RngIntElt) -> QuaternionLatticeData
{Return the quaternion lattice data for the lattice of trace zero elements in the
 Eichler order of level N in the quaternion algebra of discriminant D.}
    L, Ldual, disc_grp, to_disc, Qinv, Q, O, basis_L := get_lattice_data(D,N);
    data := New(QuaternionLatticeData);
    data`D := D;
    data`N := N;
    data`L := L;
    data`Ldual := Ldual;
    data`Q := Q;
    data`Qinv := Qinv;
    data`O := O;
    data`basis_L := basis_L;
    data`disc_grp := disc_grp;
    data`to_disc := to_disc;
    return data;
end intrinsic;

