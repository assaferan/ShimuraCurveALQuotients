// This implements the method desrcibed in [Nel15]. (Paul Nelson)
// Should also try to implement [VW] (Voight-Willis)

function W(k, x, prec)
    s := 0;
    for n in [1..prec] do
	    s +:= n^k*(n*x*KBessel(k-1, n*x) - KBessel(k, n*x));
    end for;
    return 2^(1-k)*s;
end function;

function V(k, Lbasis, t, Bound, Wprec, eps)
    RR := RealField();
    pi := Pi(RR);
    M2F := Universe(Lbasis);
    assert Type(M2F) eq AlgMat and Degree(M2F) eq 2;
    F := BaseRing(M2F);
    // We write it this way to make it easier to change the norm if we want to
    nm := map<M2F -> F | x :-> x[1,1]^2+x[2,2]^2 + x[1,2]^2 + x[2,1]^2 >;
    abs4X2 := map<M2F -> F | x :-> (x[1,1]+x[2,2])^2 + (x[1,2] - x[2,1])^2 >;
    Lgram := Matrix([[nm(x+y)-nm(x)-nm(y) : y in Lbasis] : x in Lbasis]);
    // Next line is to make sure Lgram is symmetric despite precision issues
    Lgram := 1/2*(Lgram + Transpose(Lgram));
    // Lgram := ChangeRing(Lgram, Rationals());
    L := LatticeWithGram(Lgram);
    nmBd := Bound; // if 4*|X|^2 = (a+d)^2 + (b-c)^2 < Bound, then
                   // a^2 + b^2 + c^2 + d^2 < 4*|X|^2 < Bound;
    SV := ShortVectors(L, nmBd);
    alphas := [&+[v[1][i]*Lbasis[i] : i in [1..4]] : v in SV];
    alphas_det := [alpha : alpha in alphas | Abs(Determinant(alpha) - t) lt eps];
    s := 0;
    for alpha in alphas_det do
        a,b,c,d := Explode(Eltseq(alpha));
        Xalpha := 1/2*(a+d+i*(b-c));
	    root_num := (Xalpha/ Abs(Xalpha))^k;
	    s +:= root_num * W(k, 4*pi*Abs(Xalpha), Wprec);
    end for;
    return 2*s; // because SV contains only one representative from each pair {alpha, -alpha}
end function;

function FindRdbases(D,N)
    ZZ := Integers();
    B<i,j> := QuaternionAlgebra(D);
    // need B to be indefinite
    assert exists(a){a : a in [i,j] | ZZ!(a^2) gt 0};
    if D eq 1 then
	    emb_B_BF := hom<B -> B | [1,i,j,i*j]>;
	    is_mat, _, iota := IsMatrixRing(B : Isomorphism);
    else
	    F := QuadraticField(ZZ!(a^2));
	    BF, emb_B_BF := ChangeRing(B, F);
	    is_mat, _, iota := IsMatrixRing(BF : Isomorphism);
    end if;
    assert is_mat; // F should split B
    O := Order(MaximalOrder(B), N);
    gramO := Matrix([[Norm(x+y)-Norm(x)-Norm(y) : y in Basis(O)]
		      : x in Basis(O)]);
    dual_basis := [&+[r[i]*Basis(O)[i] : i in [1..4]] : r in Rows(gramO^(-1))];
    gramO := ChangeRing(gramO, Integers());
    hnf := HermiteForm(gramO);
    basis := [&+[hnf[i,j]*dual_basis[j] : j in [1..4]] : i in [1..4]];
    Rdbases := AssociativeArray();
    for d in Divisors(D*N) do
	    Rdbases[d] := [(1/GCD(d, hnf[i,i]))*basis[i] : i in [1..4]];
	    // we actually want to embd them into M2(R)
	    Rdbases[d] := [iota(emb_B_BF(b)) : b in Rdbases[d]];
        F := BaseRing(Universe(Rdbases[d]));
        pl := RealPlaces(F)[1];
        M2R := MatrixAlgebra(RR,2);
        Rdbases[d] := [M2R | [Evaluate(x, pl) : x in Eltseq(m)] : m in Rdbases[d]];
    end for;
    return Rdbases;
end function;

// For now we only do maximal order, N = 1
function EvaluateNormAti(f, prec, Vbound, Wprec, eps, Rdbases)
    D := Level(f);
    k := Weight(f);
    s := 0;
    for n in [1..prec] do
        inner_sum := 0;
        for j->d in Divisors(D) do
            c := 1/(d*Coefficient(f,d));
            inner_sum +:= c*V(k, Rdbases[j], n/d, Vbound, Wprec, eps);
        end for;
	    s +:= Coefficient(f,n)*inner_sum;
    end for;
    return s;
end function;

// This is [Nel15, Thm 3.1]
function EvaluateNorm(D,N,f,z, prec, Vbound, Wprec, eps)
    x := Re(z);
    y := Im(z);
    sqrty := Sqrt(y);
    sigma := Matrix([[sqrty, x*sqrty^(-1)],[0,sqrty^(-1)]]);
    assert Level(f) eq D*N;
    k := Weight(f);
    Rdbases := FindRdbases(D,N);
    s := 0;
    for n in [1..prec] do
        inner_sum := 0;
        for d in Divisors(D*N) do
            c := MoebiusMu(GCD(d,N))/(d*Coefficient(f,d));
            basis := [sigma^(-1)*m*sigma : m in Rdbases[d]];
            inner_sum +:= c*V(k, basis, n/d, Vbound, Wprec, eps);
        end for;
        s +:= Coefficient(f,n)*inner_sum;
    end for;
    return y^(-k)*s;
end function;

procedure testDelta()
    delta := Basis(CuspForms(1,12))[1];
    delta_prec := 9;
    CC<i> := ComplexField();
    z := i;
    val := Evaluate(qExpansion(delta, delta_prec), Exp(2*Pi(CC)*i*z))^2;
    RdBases := FindRdbases(1,1);
    eps := 10^(-5);
    val2 := EvaluateNormAti(delta, 10, 100, 5, eps, RdBases);
    assert Abs(val - val2) lt eps;
    val3 := EvaluateNorm(1,1,delta,z,10,100,5,eps);
    assert Abs(val - val3) lt eps;
end procedure;