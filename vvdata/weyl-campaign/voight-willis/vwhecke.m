// Hecke relations, Voight-Willis eq. (17)-(18).  These are the paper's main lever: a d-dimensional
// space of forms becomes d one-dimensional eigenspaces, and a one-dimensional kernel is a far
// weaker demand on the conditioning of the system than a d-dimensional one.
//
// The paper writes (18) in the UPPER HALF-PLANE automorphy factor j(g,z) = cz+d:
//
//   K_n^h = sum_i  j(pi_i, p)^-k  j(g_i, z_{p,i})^-k  (1 - w'_{p,i})^k  (w'_{p,i})^n
//
// with the relation  sum_n K_n^h b_n = a_p f(0) = a_p b_0.  Everything else in this directory
// works in the disc, so rewrite it there.  Writing w = (z-p)/(z-pbar) for the Cayley map, so that
// 1 - w = (p - pbar)/(z - pbar), and J(g,w) for the automorphy factor of the SU(1,1) matrix
// (Phi M Phi^-1)/sqrt(det), one gets
//
//   J(g,w) = j(g,z) (1-w)/(1-w')
//
// and the three half-plane factors in (18) collapse: the (1-w) factors cancel identically and
//
//   K_n^h = sum_i  [ J(pi_i, 0) J(g_i, w_{p,i}) ]^-k  (w'_{p,i})^n .
//
// The same identity is why the Gamma rows in vw3.m are the bare  w'^n - J^k w^n  with no (1-w)^k
// weights: it is exactly eq. (11) divided by -(1-w')^k.
//
// The coset representatives are built without the group cohomology of the paper's reference.  For
// an Eichler order in an INDEFINITE algebra the one-sided class number is 1, so every left ideal of
// reduced norm p is principal; and for x of reduced norm p one has p = xbar*x, so O*x already
// contains O*p and the ideal is just O*x.  Two such x, y lie in the same coset of Gamma = O^1
// exactly when y*x^-1 lies in O, its reduced norm being 1 automatically.  So the cosets are the
// elements of reduced norm p up to left multiplication by units, and the whole problem is to
// enumerate elements of O of a given reduced norm.
//
// The reduced norm form is INDEFINITE of signature (2,2), so its "short" vectors are not a finite
// set and a coordinate box in a basis of O is the wrong shape -- it cuts the hyperboloid Nrd = p
// along directions that have nothing to do with the order.  A box of 5 gives all p+1 cosets at
// D = 6 and D = 15 for p up to 13 but only 8 of 14 at D = 22 for p = 13 and 2 of 6 at D = 26 for
// p = 5.
//
// Replace it by enumerating short vectors of the positive definite majorant FrobeniusGram of
// vw3.m, raising the bound until p+1 cosets appear.
function MajorantGram(O, mrep) return FrobeniusGram(Basis(O), mrep, Algebra(O)); end function;

// Representatives of the p+1 Hecke cosets Gamma \ {x in O : Nrd(x) = p}, for p not dividing D*N.
// Returns the representatives and the Frobenius bound (as a multiple of p) that was needed, so a
// caller can see how close the search came to its cap.
function HeckeReps(O, p : mrep := 0, prec := 60, cmax := 4096)
    A := Algebra(O);
    if mrep cmpeq 0 then mrep := FuchsianMatrixRepresentation(A : Precision := prec); end if;
    bas := Basis(O);
    L := LatticeWithGram(MajorantGram(O, mrep));
    reps := [];
    c := 8;
    while c le cmax do
        reps := [];
        for v in ShortVectors(L, SCL*c*p) do
            x := &+[ (Integers() ! v[1][i]) * bas[i] : i in [1..4] ];
            if Norm(x) ne p then continue; end if;
            if exists{ y : y in reps | (A ! y) * (A ! x)^(-1) in O } then continue; end if;
            Append(~reps, x);
            if #reps eq p+1 then return reps, c; end if;
        end for;
        c *:= 2;
    end while;
    return reps, cmax;
end function;

// The vector K_n^h of eq. (18) in the disc, together with the reduced targets, so a caller can see
// how far out of the domain the Hecke images land.
function HeckeKernelRow(G, O, DD, mrep, p0, reps, k, nterms, CC)
    Ws := SidePairingMatrices(G, DD, mrep, p0, CC);
    K := [ CC | 0 : n in [1..nterms] ];
    tgts := [];
    for x in reps do
        Px := DiscMatrix(mrep, x, CC ! p0, CC);
        jpi := DiscJ(Px, CC ! 0);
        wi := Px[1][2] / Px[2][2];
        Wg, wp := ReduceByPairing(Ws, wi, CC);
        jg := DiscJ(Wg, wi);
        Append(~tgts, <wi, wp>);
        c := (jpi*jg)^(-k);
        for n in [0..nterms-1] do K[n+1] +:= c * wp^n; end for;
    end for;
    return K, tgts;
end function;

// The eigenvalue a_p that a given expansion b satisfies, in the normalisation of (18) with the
// coset matrices scaled to determinant 1.  Used to calibrate against a known newform before the
// relation is imposed on an unknown base.
function HeckeRatio(K, bs)
    num := &+[ K[i]*bs[i] : i in [1..#K] ];
    return num / bs[1];
end function;

// One row of the linear system: sum_n K_n b_n - a_p b_0 = 0.
function HeckeRelation(K, ap)
    row := K;
    row[1] -:= ap;
    return row;
end function;
