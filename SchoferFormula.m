
intrinsic FindLambda(Q::AlgMatElt, d::RngIntElt, Order::AlgQuatOrd, basis_L::SeqEnum : bound := 10)-> BoolElt, ModTupRngElt
{Coordinate vector (in basis_L) of an optimal embedding of the discriminant-d CM order into
Order, delegating to FindLambdas (hybrid bounded-box / Embed).  Returns false if d does not embed.
The `bound' parameter is the box start bound.}
    require d gt 0: "d must be positive";
    found, lams := FindLambdas(Q, [d], Order, basis_L);
    if found then return true, lams[d]; else return false, _; end if;
end intrinsic;

intrinsic FindLambdas(Q::AlgMatElt, ds::SeqEnum[RngIntElt], Order::AlgQuatOrd, basis_L::SeqEnum : bound := 16, lambda_array := AssociativeArray())-> BoolElt, ModTupRngElt
{For each d in ds, return the integer coordinate vector v (in the basis basis_L) of an
optimal embedding of the discriminant-d CM order into Order, i.e. the element elt = sum_i
v_i/2 basis_L[i] (+ 1/2 when d = 3 mod 4) generates O cap Q(sqrt(-d)).
HYBRID strategy: PRIMARY is a BOUNDED box enumeration of the (indefinite ternary) norm form
v.Q.v = 2d -- deterministic and, unlike Magma's Embed, it never hangs; the box bound starts at
`bound' and doubles only up to a hard cap, so it can never blow up (the old unbounded box started at
max_d/2 and 80 GB-exploded).  FALLBACK is Magma's Embed, used only for discs whose optimal vector
exceeds the box cap (large quadratic-CM discs).  This keeps the small rational discs -- where Embed is
the one that intermittently hangs on unlucky ShimuraCurveLattice order representatives -- on the safe
box path.  Returns true with the lambdas iff all d found.}
    require &and[d gt 0 : d in ds]: "All ds must be positive";
    lambdas := lambda_array;
    QZ := ChangeRing(Q, Integers());
    n := Nrows(QZ);
    A := Algebra(Order);
    // express algebra elements in the (trace-zero) basis basis_L
    ML := Matrix([Eltseq(A!b) : b in basis_L]);
    coordsL := func< el | Solution(ML, Vector(Eltseq(A!el))) >;
    R<x> := PolynomialRing(Rationals());
    B_O := Basis(Order);
    Mat_O := Matrix([Eltseq(B_O[i]) : i in [1..#B_O]]);
    Sprime := HorizontalJoin(IdentityMatrix(Integers(),2), ZeroMatrix(Integers(),2));
    e1col := Eltseq(Solution(Mat_O, Vector(Rationals(),[1,0,0,0])));

    // does the integer coordinate vector v give an OPTIMAL embedding of the disc-(-d) order?
    // elt = sum v_i/2 basis_L[i] (+1/2 for d = 3 mod 4) must generate O cap Q(sqrt(-d)); the
    // Smith-form test below checks that {1, elt} is part of a Z-basis of O.
    is_optimal := function(v, d)
        if d mod 4 eq 3 then
            elt := &+[v[i]/2*basis_L[i] : i in [1..#basis_L]] + Parent(basis_L[1])!1/2;
        else
            elt := &+[v[i]/2*basis_L[i] : i in [1..#basis_L]];
        end if;
        CoB := Matrix([e1col, Eltseq(Solution(Mat_O, Vector(Rationals(), Eltseq(elt))))]);
        CoBZZ := ChangeRing(Denominator(CoB)*CoB, Integers());
        return SmithForm(CoBZZ) eq Sprime;
    end function;

    // ---- PRIMARY: bounded box enumeration of v.Q.v = 2d (deterministic; never hangs) ----
    // Trace-zero space is 3-dimensional for a quaternion order; evaluate the ternary form directly.
    if n eq 3 then
        twice := {Integers()| 2*d : d in ds};
        q11 := QZ[1,1]; q22 := QZ[2,2]; q33 := QZ[3,3];
        q12 := QZ[1,2]+QZ[2,1]; q13 := QZ[1,3]+QZ[3,1]; q23 := QZ[2,3]+QZ[3,2];
        // cap must clear the LARGEST box bound any required d needs. This depends on the
        // (non-canonical) order representative Magma hands us, so a borderline d found right
        // at the cap on a lucky run can spill past it on an unlucky one -- the box then misses
        // it and the flaky Embed fallback may fail to recover it (observed intermittently for
        // d = 280 in X0(10,19), found at bound 128 on lucky runs). Keep the cap generous so the
        // deterministic box always wins; the loop exits as soon as all d are found, so the
        // higher cap costs nothing on the common (small-bound) case.
        cap := 512;
        bd := Maximum(bound, 8);
        while (bd le cap) and not (Set(ds) subset Keys(lambdas)) do
            for a in [-bd..bd] do
                for b in [-bd..bd] do
                    pab := q11*a^2 + q22*b^2 + q12*a*b;
                    qc  := q13*a + q23*b;
                    for c in [-bd..bd] do
                        nv := pab + q33*c^2 + qc*c;
                        if nv in twice then
                            d := nv div 2;
                            if d notin Keys(lambdas) then
                                v := Vector([Integers()| a, b, c]);
                                if is_optimal(v, d) then
                                    lambdas[d] := v;
                                    vprintf ShimuraQuotients, 3 : "\t  box: d = %o at bound <= %o\n", d, bd;
                                end if;
                            end if;
                        end if;
                    end for;
                end for;
            end for;
            bd *:= 2;
        end while;
    end if;

    // ---- FALLBACK: Magma's Embed for any d the bounded box did not reach ----
    //   d = 3 mod 4: order Z[(1+sqrt(-d))/2] (disc -d), generator (1+sqrt(-d))/2 -> lam = 2*mu-1.
    //   d = 0 mod 4: order Z[sqrt(-d/4)]    (disc -d), generator sqrt(-d/4)     -> lam = 2*mu.
    // (Embedding Z[sqrt(-d)] (disc -4d) instead would make elt = lam/2 fall outside Order.)
    for d in ds do
        if d in Keys(lambdas) then continue; end if;
        vprintf ShimuraQuotients, 2 : "\t  box missed d = %o; Embed fallback (may hang -- see FindLambdas)\n", d;
        embeds := true; lam := A!0;
        try
            if d mod 4 eq 3 then
                lam := 2*(A!Embed(MaximalOrder(NumberField(x^2+d)), Order)) - 1;
            else
                lam := 2*(A!Embed(EquationOrder(NumberField(x^2+(d div 4))), Order));
            end if;
        catch e
            embeds := false;                                // K does not embed -> no lambda for this d
        end try;
        if not embeds then continue; end if;
        cv := coordsL(lam);
        if not &and[IsCoercible(Integers(), cc) : cc in Eltseq(cv)] then continue; end if;
        v := Vector([Integers()| cc : cc in Eltseq(cv)]);
        if (v*QZ, v) ne 2*d then continue; end if;          // self-check: correct norm
        if is_optimal(v, d) then lambdas[d] := v; end if;
    end for;

    if Set(ds) eq Keys(lambdas) then
        return true, lambdas;
    else
        vprint ShimuraQuotients, 2 : "Could not find all lambdas";
        vprint ShimuraQuotients, 2 : Set(ds) diff Keys(lambdas);
        return false, lambdas; //return the partial progress
    end if;
end intrinsic;

intrinsic ElementOfNorm(Q::AlgMatElt, d::RngIntElt, Order::AlgQuatOrd, basis_L::SeqEnum) -> ModTupRngElt
{Coordinate vector (in basis_L) of an optimal embedding of the disc-d CM order into Order,
via FindLambda (Magma's Embed).  Errors if d does not embed.}
    require d gt 0: "d must be a positive norm";
    found_lambda, lambda := FindLambda(Q, d, Order, basis_L);
    if not found_lambda then
        error Sprintf("ElementOfNorm: no optimal embedding of discriminant %o into the order", d);
    end if;
    return lambda;
end intrinsic;

intrinsic ElementsOfNorm(Q::AlgMatElt, ds::SeqEnum[RngIntElt], Order::AlgQuatOrd, basis_L::SeqEnum) -> ModTupRngElt
{For each d in ds, return the coordinate vector (in basis_L) of an optimal embedding of the
disc-d CM order into Order, via FindLambdas (Magma's Embed).  Errors if some d does not embed.}
    require &and[d gt 0 : d in ds]: "All ds must be positive";
    vprintf ShimuraQuotients, 2 : "\n\tFinding lambdas for norms in %o (optimal embeddings)...", ds;
    found_lambdas, lambdas := FindLambdas(Q, ds, Order, basis_L);
    if not found_lambdas then
        error Sprintf("ElementsOfNorm: no optimal embedding found for discriminant(s) %o",
                      Set(ds) diff Keys(lambdas));
    end if;
    vprintf ShimuraQuotients, 2 : "Found lambdas.";
    return lambdas;
end intrinsic;

intrinsic VerticalJoinList(mats::List)->.
    {}
    m := mats[1];
    for i in [2..#mats] do
        m := VerticalJoin(m, mats[i]);
    end for;
    return m;
end intrinsic;

function my_legendre_symbol(alpha, p)
    return LegendreSymbol(Integers()!(GF(p)!alpha),p);
end function;

// W_{m,p}
// L should be Lminus
function Wpoly(m,p,mu,L,K,Q)
    _<sqrtp> := K;
    F := QNF();
    S := BasisMatrix(L)*Q*Transpose(BasisMatrix(L));
    n := Rank(L);
    Lnf := NumberFieldLatticeWithGram(ChangeRing(S,F));
    // l is the sequence of exponents
    assert p ne 2; // take care of p = 2 later
    bases, Jblocks, exps := JordanDecomposition(Lnf,p*Integers(F)); 
    l := &cat[[e : j in [1..Nrows(Jblocks[i])]] : i->e in exps];
    eps := &cat[[Rationals() | 1/2 * x / p^exps[i] : x in Diagonal(b)]  : i->b in Jblocks]; // so that S is equivalent to (2 eps_1 p^{l_1},..., 2 eps_n p^{l_n})
    assert &and[Valuation(e,p) eq 0 : e in eps];
    B := ChangeRing(VerticalJoinList(bases), Rationals());
    mu_wrt_L := Solution(ChangeRing(BasisMatrix(L),Rationals()), ChangeRing(mu, Rationals()));
    Q_mu := 1/2*(mu_wrt_L * ChangeRing(S,Rationals()), mu_wrt_L);
    R<x> := PolynomialRing(K);
    if not IsIntegral(m - Q_mu) then
        return R!0;
    end if;
    mu_wrt_B := mu_wrt_L*B^(-1);
    H_mu := {i : i in [1..n] | Valuation(mu_wrt_B[i],p) ge 0};
    vals := [l[i] + Valuation(mu_wrt_B[i], p) : i in [1..n] | i notin H_mu];
    if IsEmpty(vals) then
        K_0 := Infinity();
    else
        K_0 := Minimum(vals);
    end if;
    L_mu := func<k | {i : i in H_mu | IsOdd(l[i] - k) and (l[i] - k lt 0)}>;
    l_mu := func<k | #L_mu(k)>;
    // we compute twice d_mu for technical reasons
    d2_mu := func<k | 2*k + &+[Minimum(l[i]-k, 0) : i in H_mu]>;
    eps_mu := func<k | LegendreSymbol(-1,p)^(l_mu(k) div 2) * &*[Integers() | my_legendre_symbol(eps[i],p) : i in L_mu(k)]>;
    f_1 := function(x)
        a, alpha := Valuation(x,p);
        return IsEven(l_mu(a+1)) select -1/p else my_legendre_symbol(alpha,p) / sqrtp;
    end function;
    t_mu := m - &+[Rationals() | eps[i]*p^l[i]*mu_wrt_B[i]^2 : i in [1..n] | i notin H_mu];
    a := Valuation(t_mu, p);
    
    if a lt K_0 then
        ret := 1;
        ret +:= (1 - 1/p)*&+[R | eps_mu(k)*sqrtp^d2_mu(k)*x^k : k in [1..a] | IsEven(l_mu(k))];
        ret +:= eps_mu(a+1)*f_1(t_mu)*sqrtp^d2_mu(a+1)*x^(a+1);
    else
        ret := 1;
        ret +:= (1 - 1/p)*&+[R | eps_mu(k)*sqrtp^d2_mu(k)*x^k : k in [1..K_0] | IsEven(l_mu(k))];
    end if;
    return ret;
end function;

// This is for p = 2
function Wpoly2(m,mu,L,K,Q)
    p := 2;
    Zp := pAdicRing(p);
    _<sqrtp> := K;
    F := QNF();
    S := BasisMatrix(L)*Q*Transpose(BasisMatrix(L));
    n := Rank(L);
    Lnf := NumberFieldLatticeWithGram(ChangeRing(S,F));
    // l is the sequence of exponents
    assert p eq 2;
    bases, Jblocks, exps := JordanDecomposition(Lnf,p*Integers(F)); 
    bases := [* ChangeRing(B, Rationals()) : B in bases *];
    Jblocks := [* ChangeRing(J, Rationals()) : J in Jblocks *];
    bases := [* ChangeRing(B, Zp) : B in bases *];
    Jblocks := [* ChangeRing(J, Zp) : J in Jblocks *];
    l_list := [];
    m_list := [];
    n_list := [];
    eps := [];
    mu_indices := [];
    // For these we record the first index of two, so that (i, i+1) are the indices
    mu_prime_indices := [];
    mu_prime_prime_indices := [];
    // It seems like eps_prime and eps_prime_prime can always be taken to be 1
    // eps_prime := [];
    // eps_prime_prime := [];
    row_ind := 0;
    for i->Jblock in Jblocks do
        if Nrows(Jblock) eq 1 then
            row_ind +:= 1;
            Append(~l_list, exps[i]);
            Append(~eps, Jblock[1,1] / p^exps[i]);
            Append(~mu_indices, row_ind);
        end if;
        for j in [2..Nrows(Jblock)] do
            row_ind +:= 1;
            b := Jblock[j-1,j] / p^(exps[i]);
            if IsWeaklyZero(b) then
                Append(~l_list, exps[i]);
                Append(~eps, Jblock[j-1,j-1] / p^exps[i]);
                Append(~mu_indices, row_ind);
                if (j eq Nrows(Jblock)) then
                    row_ind +:= 1;
                    Append(~l_list, exps[i]);
                    Append(~eps, Jblock[j,j] / p^exps[i]);
                    Append(~mu_indices, row_ind);
                end if;
                continue; 
            end if;
            a := Jblock[j-1,j-1] / p^(exps[i]);
            d := Jblock[j,j] / p^(exps[i]);
            disc := b^2 - a*d;
            if (Integers(8)!(Integers()!disc) eq 5) then 
                disc +:= 2*a; 
                Append(~n_list, exps[i]);
                Append(~mu_prime_prime_indices, row_ind);
                aniso := true;
            else   
                Append(~m_list, exps[i]);
                Append(~mu_prime_indices, row_ind);
                aniso := false;
            end if;
            is_sqr, sqrt_disc := IsSquare(disc);
            assert is_sqr;
            if IsWeaklyZero(a) then
                B := Matrix(Zp, [[1,0],[0,1]]);
            else
                // solving the quadratic
                x1 := (-b + sqrt_disc) div a;
                x2 := (-b - sqrt_disc) div a;
                if aniso then
                    inner_product := 2-(2*disc div a);
                    x := Zp!inner_product;
                    beta := Sqrt((4-x^2)/3);
                    alpha := (-beta + x)/2;
                    // if v_1, v2 are a basis for [[2,1],[1,2]], then v_1, alpha v1 + beta v2 are a basis for [[2,x],[x,2]]
                    get_to_x := Matrix([[1,0],[alpha, beta]]);
                    assert &and[IsWeaklyZero(e) : e in Eltseq(get_to_x * Matrix([[2,1],[1,2]]) * Transpose(get_to_x) - Matrix([[2,x],[x,2]]))];
                    B := Matrix([[x1, 1], [x2, 1]]);
                    B := get_to_x^(-1)*B;
                else
                    z2 := (-a)/(2*disc); // constant to get scalar product equal to 1 
                    // Change of basis matrix
                    B := Matrix([[x1, 1], [z2*x2, z2]]);
                end if;
            end if;
            cans := [SymmetricMatrix([0,1,0]), SymmetricMatrix([2,1,2])];
            assert B*Matrix([[a,b],[b,d]])*Transpose(B) in cans;
            // !!!! This stopped working in V2.29 !!!!
            // B_big := Parent(Jblock)!1;
            B_big := IdentityMatrix(BaseRing(Jblock),Nrows(Jblock));
            B_big[j-1,j-1] := B[1,1];
            B_big[j-1,j] := B[1,2];
            B_big[j,j-1] := B[2,1];
            B_big[j,j] := B[2,2];
            bases[i] := B_big * bases[i];
            Jblocks[i] := B_big * Jblocks[i] * Transpose(B_big);
            row_ind +:= 1;
        end for;
    end for;
    eps_prime := [1 : m in m_list];
    eps_prime_prime := [1 : n in n_list];
    H := #l_list;
    M := #m_list;
    N := #n_list;
    assert n eq H + 2*M + 2*N;
    assert &and[Valuation(e) eq 0 : e in eps];
    B := ChangeRing(VerticalJoinList(bases), Rationals());
    mu_wrt_L := Solution(ChangeRing(BasisMatrix(L),Rationals()), ChangeRing(mu, Rationals()));
    Q_mu := 1/2*(mu_wrt_L * ChangeRing(S,Rationals()), mu_wrt_L);
    R<x> := PolynomialRing(K);
    if not IsIntegral(m - Q_mu) then
        return R!0;
    end if;
    mu_wrt_B := mu_wrt_L*B^(-1);
    mu_list := [mu_wrt_B[i] : i in mu_indices];
    mu_prime_list := [[mu_wrt_B[i], mu_wrt_B[i+1]] : i in mu_prime_indices];
    mu_prime_prime_list := [[mu_wrt_B[i], mu_wrt_B[i+1]] : i in mu_prime_prime_indices];
    M_mu := {i : i in [1..M] | Valuation(mu_prime_list[i][1],p) ge 0 and Valuation(mu_prime_list[i][2],p) ge 0};
    N_mu := {i : i in [1..N] | Valuation(mu_prime_prime_list[i][1],p) ge 0 and Valuation(mu_prime_prime_list[i][2],p) ge 0};
    H_mu := {i : i in [1..H] | Valuation(mu_list[i],p) ge 0};

    L_mu := func<k | {i : i in H_mu | IsOdd(l_list[i] - k) and (l_list[i] - k lt 0)}>;
    l_mu := func<k | #L_mu(k)>;
    // we compute twice d_mu for technical reasons
    d2_mu := func<k | 2*k + &+[Integers()|Minimum(l_list[i]-k, 0) : i in H_mu] + 2*&+[Integers()|Minimum(m_list[i]-k,0) : i in M_mu] + 2*&+[Integers()|Minimum(n_list[i]-k,0) : i in N_mu]>;

    vals := [l_list[i] + Valuation(mu_list[i], p) : i in [1..H] | i notin H_mu and Valuation(mu_list[i],p) lt -1];
    vals cat:= [l_list[i] : i in [1..H] | i notin H_mu and Valuation(mu_list[i],p) eq -1];
    vals cat:= [m_list[i] + Minimum(Valuation(mu_prime_list[i][1], p), Valuation(mu_prime_list[i][2], p)) : i in [1..M] | i notin M_mu];
    vals cat:= [n_list[i] + Minimum(Valuation(mu_prime_prime_list[i][1], p), Valuation(mu_prime_prime_list[i][2], p)) : i in [1..N] | i notin N_mu];
    if IsEmpty(vals) then
        K_0 := Infinity();
    else
        K_0 := Minimum(vals);
    end if;
    
    p_mu := func<k | (-1)^(&+[Integers()|Minimum(n_list[j] - k, 0) : j in N_mu])>;
    eps_mu := func<k | &*[Zp|eps[h] : h in L_mu(k)]>;

    delta_mu := func<k | exists(h){h : h in H_mu | l_list[h] eq k} select 0 else 1>;

    two_block := func< x | x[1]^2 + x[1]*x[2] + x[2]^2>;

    Q_prime_mu := &+[Rationals() | eps[i]*p^(l_list[i]-1)*mu_list[i]^2 : i in [1..H] | i notin H_mu];
    Q_prime_mu +:= &+[Rationals() | eps_prime[i]*p^m_list[i]*(&* mu_prime_list[i]) : i in [1..M] | i notin M_mu];
    Q_prime_mu +:= &+[Rationals() | eps_prime_prime[i]*p^n_list[i]*two_block(mu_prime_prime_list[i]) : i in [1..N] | i notin N_mu];

    t_mu := m - Q_prime_mu;
    a := Valuation(t_mu, p);

    nu := func< k | t_mu*p^(3-k) - &+[Zp|eps[h] : h in H_mu | l_list[h] lt k]>;

    psi_char := func<k | (Valuation(nu(k)) ge 2) select (-1)^(Integers()!(GF(2))!(nu(k)/ 4)) else 0>;

    K_0 := Minimum(K_0, a+3);

    
    ret := 1;
    ret +:= &+[R | delta_mu(k)*p_mu(k)*sqrtp^(d2_mu(k)-3)*KroneckerSymbol(2,Integers()!(eps_mu(k)*nu(k)))*x^k : k in [1..K_0] | IsOdd(l_mu(k))];
    ret +:= &+[R | delta_mu(k)*p_mu(k)*sqrtp^(d2_mu(k)-2)*KroneckerSymbol(2,Integers()!(eps_mu(k)))*psi_char(k)*x^k : k in [1..K_0] | IsEven(l_mu(k))];

    return ret;
end function;

function Wpoly_scaled(m,p,mu,L,Q : scaled := true)
    S := BasisMatrix(L)*Q*Transpose(BasisMatrix(L));
    D := Determinant(S);
    vpD := Valuation(D,p);
    K<sqrtp> := QuadraticField(p);
    scale := sqrtp^(-vpD);
    euler := p eq 2 select Wpoly2(m,mu,L,K,Q) else Wpoly(m,p,mu,L,K,Q);
    assert CanChangeRing(euler, Rationals());
    return (scaled) select scale*euler else ChangeRing(euler, Rationals()), p^(-vpD);
end function;

function W(m,p,mu,L,Q)
    Wpoly := Wpoly_scaled(m,p,mu,L,Q);
    _<sqrtp> := BaseRing(Wpoly);
    n := Rank(L);
    // s0 := n/2 - 1;
    // s := -s0;
    s2 := 2 - n; // s2 = 2*s always integral
    return Evaluate(Wpoly, sqrtp^(-s2));
end function;

function get_Wpolys(m,mu,Lminus,Q, Sm_mu : scaled := true)
    if scaled then 
        Wpolys := [* Wpoly_scaled(m,p,mu,Lminus,Q) : p in Sm_mu *];
        wpolyseval := [* Evaluate(Wpolys[i],1) : i in [1..#Wpolys] *];
        return Wpolys, wpolyseval;
    end if;

    Wpolys := [];
    scales_sqr := [Rationals() | ];
    wpolyseval := [];
    for p in Sm_mu do
        Wpol, scale_sqr :=  Wpoly_scaled(m,p,mu,Lminus,Q : scaled := false);
        Append(~Wpolys, Wpol);
        Append(~scales_sqr, scale_sqr);
        Append(~wpolyseval, Evaluate(Wpol, 1));
    end for;
    
    return Wpolys, wpolyseval, scales_sqr;
end function;

function find_v(F, places, sqrtps)
    // assert exists(v){v : v in places | &and[Evaluate(F!sqrtp, v) gt 0 : sqrtp in sqrtps]};
    pos := false;
    for v in places do
        pos := true;
        for sqrtp in sqrtps do
            if Evaluate(F!sqrtp, v) lt 0 then pos := false; break; end if;
        end for;
        if pos then v_pos := v; break; end if;
    end for;
    assert pos;
    return v_pos;
end function;

function get_kappa_minus_squared(d, Wpolys, Wpol, Sm_mu, i, scales_sqr)
    W_prod := &*[ Rationals() | Evaluate(Wpolys[j],1) : j in [1..#Sm_mu] | j ne i];
    W_prod *:= Evaluate(Derivative(Wpol),1); // this should be multiplied by -log(p_prime)

    kron_prod := &*[Rationals() | 1 - Evaluate(KroneckerCharacter(d),p)/p : p in Sm_mu];

    h := ClassNumber(d);
    w := #UnitGroup(QuadraticField(d));

    scale_sqr := &*scales_sqr;
    
    W_kron := W_prod / kron_prod;
    // km_sqr := -d*scale_sqr*(w*W_kron / h)^2;
    // Using Yang's code to try to work the non-maximal case
    // !!! Not sure why this works !!!
    _, f := SquarefreeFactorization(d div FundamentalDiscriminant(d));
    km_sqr := -d*scale_sqr*(w*W_kron / h)^2 / f^2;
    km_sign := -Sign(W_kron);

    return km_sqr, km_sign;
   
end function;

// returns x,y such that the answer is x logy
function kappaminus(mu, m, Lminus, Q, d)
    // error if m eq 0, "Not implemented for m eq 0 at CM point!\n", d;  
    if m eq 0 then
        printf "\nWarning: Adding 0 for m eq 0 at CM point!\n";
        return 0, 1;
    end if;
    Bminus := BasisMatrix(Lminus);
    Delta := Determinant(Bminus*Q*Transpose(Bminus));
    
    Sm_mu := {p : p in PrimeDivisors(Delta)} join {p : p in PrimeDivisors(Numerator(m))};
    Sm_mu := [p : p in Sm_mu];
    
    vprintf ShimuraQuotients, 5: "\t\t\tSm_mu = %o\n", Sm_mu;

    Wpolys, wpolyseval, scales_sqr := get_Wpolys(m,mu,Lminus,Q, Sm_mu : scaled := false);
   
    assert exists(i){i : i in [1..#Sm_mu] | wpolyseval[i] eq 0};
    p_prime := Sm_mu[i];
    if exists(j){j : j in [1..#Sm_mu] | wpolyseval[j] eq 0 and j ne i} then
        return 0, p_prime;
    end if;
    Wpol := Wpolys[i];

    // F, sqrtps := get_field(Wpolys);

    // ret := get_kappa_minus(F, d, Wpolys, Wpol, Sm_mu, i, sqrtps);
    // ret_squared, ret_sign := get_kappa_minus_squared(F, d, Wpolys, Wpol, Sm_mu, i, sqrtps);
    ret_squared, ret_sign := get_kappa_minus_squared(d, Wpolys, Wpol, Sm_mu, i, scales_sqr);
    // assert ret_squared eq ret^2 and ret_sign eq Sign(ret);
    is_sqr, ret := IsSquare(ret_squared);
    assert is_sqr;
    ret := ret_sign*AbsoluteValue(ret);

    vprintf ShimuraQuotients, 5 : "\t\t\tadding %o log %o\n", -ret, p_prime;
    return -ret, p_prime; // to get x logy instead of -xlogy
    // return p_prime^(-ret);
end function;

function kappaminuszero(D,N,d)
    log_coeffs := AssociativeArray();
    for p in PrimeDivisors(D div GCD(d,D)) do
        log_coeffs[p] := (p-1)/(p+1);
    end for;
    for p in PrimeDivisors(N div GCD(d,N)) do
        log_coeffs[p] := 1;
    end for;
    // !! TODO := think about precision
    RR := RealField();
    pi := Pi(RR);
    gamma :=  EulerGamma(RR);
    chi := KroneckerCharacter(d);
    mu := #UnitGroup(QuadraticField(d));
    h := ClassNumber(d);
    chowla_selberg := &+[chi(a)*mu*Log(Gamma(a/AbsoluteValue(d)))/h : a in [1..AbsoluteValue(d)-1]];
    chowla_selberg +:= Log(4*pi) - 3*Log(AbsoluteValue(d)) + gamma;
    return log_coeffs, chowla_selberg;
end function;

// Computes kappa0(m) in Schofer's formula
intrinsic Kappa0(m::RngIntElt, d::RngIntElt, Q::AlgMatElt, lambda_v::ModTupRngElt) -> LogSm
{Computing coefficients Kappa0(m) in Schofers formula}
    return Kappa(Parent(lambda_v)!0,Rationals()!m,d,Q,lambda_v);
end intrinsic;

intrinsic Kappa(gamma::ModTupRngElt, m::FldRatElt, d::RngIntElt, Q::AlgMatElt, lambda_v::ModTupRngElt) -> LogSm
{Computing coefficients Kappa(gamma, m) in Schofers formula.}
    vprintf ShimuraQuotients, 4: "\n\t\t\tKappa_%o of %o", gamma, m;
    Qrat := ChangeRing(Q, Rationals());
    Q := ChangeRing(Q, Integers());
    
    c_Lplus := Content(lambda_v);
    Lplus := RSpaceWithBasis(Matrix(lambda_v div c_Lplus));
    Lminus := Kernel(Transpose(Matrix(lambda_v*Q)));
    L := RSpaceWithBasis(IdentityMatrix(Integers(),3));
    L_quo, L_quo_map := L / (Lplus + Lminus);

    lambda_rat := ChangeRing(lambda_v, Rationals());
    gamma_rat := ChangeRing(gamma, Rationals());
    c_gamma_plus := ((gamma_rat*Qrat, lambda_rat)/(lambda_rat*Qrat,lambda_rat));
    gamma_plus:= c_gamma_plus*lambda_rat;
    gamma_minus := gamma_rat - gamma_plus;
    log_coeffs := LogSum();
    // This is the condition from Yang code, if we have a vector with Q(x) = m
    Yang_tt := false;
    for mu_bar in L_quo do
        mu := mu_bar@@L_quo_map;
        c_mu_plus := ((mu*Q, lambda_v)/(lambda_v*Q,lambda_v));
        mu_plus:= c_mu_plus*ChangeRing(lambda_v, Rationals());
        mu_minus := mu - mu_plus;
        // finding the possible range of x in mu_plus + L_plus
        // use that mu_plus = c_mu_plus * lambda, L_plus = Z * c_Lplus^(-1) * lambda
        // gamma_plus = c_gamma_plus*lambda
        // that <lambda,lambda> = -2d, and we only need x with <x,x> <= 2m
        // so if x = c_gamma_plus + c_mu_plus + c_Lplus^(-1)*k, we need only those with
        // (c_gamma_plus + c_mu_plus + c_Lplus^(-1)*k)^2 le m/(-d)
        // thus k is between the following bounds
        sqr_bd := m/(-d);
        
        lb := Ceiling((-Sqrt(sqr_bd) - c_mu_plus - c_gamma_plus)*c_Lplus);
        ub := Floor((Sqrt(sqr_bd) - c_mu_plus - c_gamma_plus)*c_Lplus);
        for k in [lb..ub] do
            x := (c_gamma_plus + c_mu_plus + k * c_Lplus^(-1)) * lambda_rat;
            assert (m - (x*Qrat,x)/2) ge 0;
            // if (m - (x*Qrat,x)/2) lt 0 then printf "skipping...\n"; continue; end if;
            vprintf ShimuraQuotients, 5: "\n\t\t\tmu_minus = %o, m - Q(x) = %o\n", gamma_minus + mu_minus, m - (x*ChangeRing(Q,Rationals()),x)/2;
            norm_mu_minus := ((gamma_minus + mu_minus)*Qrat, gamma_minus + mu_minus)/2;
            vprintf ShimuraQuotients, 5: "\t\t\tQ(mu_minus) = %o, Q(mu_minus) - m + Q(x) = %o\n", norm_mu_minus, norm_mu_minus - m + (x*ChangeRing(Q,Rationals()),x)/2;
            if (m - (x*Qrat,x)/2 eq 0) then // and (gamma_minus + mu_minus ne 0) then // This condition is for Chowla-Selberg constant
                if (gamma ne 0) then
                    Yang_tt := true;
                else
                    m0, m_cond := SquareFreeFactorization(Integers()!m);
                    fac := Factorization(m_cond);
                    for pe in fac do
                        p,e := Explode(pe);
                        log_coeffs -:= LogSum(Rationals()!2*e,p);
                    end for;
                end if;
            else
                a, p := kappaminus(gamma_minus + mu_minus, m - (x*Qrat,x)/2, Lminus, Q, d);
                log_coeffs +:= LogSum(Rationals()!a,p);
            end if;
        end for;
    end for;

    // trying to imitate Yang's code
    // !!! Don't know why this is working !!!
    if Yang_tt then
        d0 := FundamentalDiscriminant(d);
        f2 := d div d0;
        is_pp, p, e2 := IsPrimePower(f2);
        if is_pp then
            e := e2 div 2;
            log_coeffs +:= LogSum(-4*p^(1-e)/(p-KroneckerSymbol(d0,p)),p);
        end if;
    end if;

    return log_coeffs;
end intrinsic;

intrinsic SchoferFormula(fs::SeqEnum[RngSerLaurElt], d::RngIntElt, Q::AlgMatElt, lambda::ModTupRngElt, scale::FldRatElt) -> SeqEnum[LogSm]
{Assuming that fs are q-expansions of oo-weakly holomorphic modular forms at oo, 
 returns the log of the absolute value of Psi_F_f for every f in fs at the CM point with CM d.
 Here Q is the Gram matrix of the lattice L and lambda is a vecotr of norm -d.}
    ns := [-Valuation(f) : f in fs];
    n := Maximum(ns);
    log_coeffs := [LogSum() : f in fs];
    for m in [1..n] do
        if &and[Coefficient(f, -m) eq 0 : f in fs] then continue; end if;
        log_coeffs_m := Kappa0(m,d,Q,lambda);
        vprintf ShimuraQuotients, 5 : "\t\t";
        vprintf ShimuraQuotients, 4 : " is %o", log_coeffs_m;
        for i->f in fs do
            log_coeffs[i] +:= Coefficient(f,-m)*log_coeffs_m;
        end for;
    end for;

    // rescaling

    for i in [1..#fs] do
        log_coeffs[i] := scale * log_coeffs[i];
    end for;

    return log_coeffs;
end intrinsic;

function SchoferFormula0(fs_0, d, Q, lambda_v, scale, M, disc_grp, to_disc)

    log_coeffs := [LogSum() : f in fs_0];

    ns := [-Valuation(f) : f in fs_0];
    n := Maximum(ns);
    
    // computing norms of elements in the discriminant group
    mod_M_to_vecs := AssociativeArray([0..M-1]);
    for j in [0..M-1] do
        mod_M_to_vecs[j] := [];
    end for;
    for eta in disc_grp do
        v := ChangeRing(eta@@to_disc,Rationals());
        norm_v := (v*Q,v)/(2*M);
        if not IsIntegral(norm_v) then continue; end if;
        norm_mod_M := Integers()!norm_v mod M;
        Append(~mod_M_to_vecs[norm_mod_M], eta);
    end for;

    for mM in [1..n] do
        if &and[Coefficient(f, -mM) eq 0 : f in fs_0] then continue; end if;
        gammas:= [1/M*ChangeRing(gammaM@@to_disc, Rationals()) : gammaM in mod_M_to_vecs[mM mod M]];
        log_coeffs_m := &+([Kappa(gamma,mM/M,d,Q,lambda_v) : gamma in gammas] cat [LogSum()]);
        vprintf ShimuraQuotients, 4 : " is %o", log_coeffs_m;
        for i->f in fs_0 do
            log_coeffs[i] +:= Coefficient(f,-mM)*log_coeffs_m;
        end for;
    end for;

    // rescaling

    for i in [1..#fs_0] do
        log_coeffs[i] := scale * log_coeffs[i];
    end for;

    return log_coeffs;
end function;

intrinsic SchoferFormula(f::RngSerLaurElt, d::RngIntElt, Q::AlgMatElt, lambda::ModTupRngElt, scale::FldRatElt) -> LogSm
{Assuming that f is the q-expansions of a oo-weakly holomorphic modular form at oo, 
 returns the log of the absolute value of Psi_F_f at the CM point with CM d.
 Here Q is the Gram matrix of the lattice L and lambda is a vecotr of norm -d.}
    return SchoferFormula([f], d, Q, lambda, scale)[1];
end intrinsic;

intrinsic ScaleForSchofer(d::RngIntElt, D::RngIntElt, N::RngIntElt) -> FldRatElt
{Return the scaling factor in Schofer formula for CM(d) on X*(D,N).}
    D0 := D div 2^Valuation(D,2);
    M := 4*D0;
    
    OK := MaximalOrder(QuadraticField(d));
    is_sqr, cond := IsSquare(d div Discriminant(OK));
    assert is_sqr;
    // require cond eq 1 : "Not implemented for non-maximal orders!";
    O := sub<OK | cond>;
    n_d := NumberOfOptimalEmbeddings(O, D, N);
    require n_d gt 0 : "Curve does not have a CM point of discirminant d!";
    W_size := 2^#PrimeDivisors(D*N);

    // Not sure?? Think this what happens to the number of CM points on the full quotient
    /*
    sqfree, sq := SquarefreeFactorization(d);
    Ogg_condition := (cond eq 1) or (((sqfree mod 4 eq 1) and (cond eq 2)));
    if ((D*N) mod sqfree eq 0) and Ogg_condition then
        W_size div:= 2;
    end if;
    */
    // This follows from Ogg's description of the fixed points 
    // of Atkin-Lehner w_m
    Ogg_condition := ((d eq -4) and IsEven(D*N)) or
                     ((d mod 4 eq 0) and ((D*N mod (d div 4)) eq 0)) or
                     ((d mod 4 eq 1) and (D*N mod d eq 0));
    if Ogg_condition then
        W_size div:= 2;
    end if;

    scale := -n_d / (4*W_size);

    return scale;
end intrinsic;

intrinsic SchoferFormula(f::RngSerLaurElt, d::RngIntElt, D::RngIntElt, N::RngIntElt, Ldata::QuaternionLatticeData : Lambda := false) -> LogSm
{Assuming that f is the q-expansions of a oo-weakly holomorphic modular form at oo, 
 returns the log of the absolute value of Psi_F_f at the CM point with CM d.}
    // _,_,_,_,_,Q,O,basis_L := ShimuraCurveLattice(D,N);
    Q := Ldata`Q;
    O := Ldata`O;
    basis_L := Ldata`basis_L;

    scale := ScaleForSchofer(d,D,N);

     if Type(Lambda) eq BoolElt then 
        lambda := ElementOfNorm(Q,-d, O, basis_L);
    else
        lambda := Lambda;
    end if;

    return SchoferFormula(f, d, Q, lambda, scale);
end intrinsic;


// !! TODO - cache the Kappa0's or do it for a bunch of fs simultaneously
// We use a combination of the two versions of Schofer's formula from [GY] and [Err]
// We write sum log|psi|^2 = -|CM(d)|/4 * sum c_m kappa(-m)
// Note that in [GY] there is no square on the lhs, and 
// in [Err] there is no division by 4 on the rhs,
// but this seems to match with the examples in [Err] !?
intrinsic SchoferFormula(etas::SeqEnum[EtaQuot], d::RngIntElt, D::RngIntElt, N::RngIntElt, Ldata::QuaternionLatticeData : Lambda := false) -> SeqEnum[LogSm]
{Return the log of the absolute value of Psi_F_f for every f in fs at the CM point with CM d.}
    // _,_,disc_grp,to_disc,_, Q, O, basis_L := ShimuraCurveLattice(D,N);
    Q := Ldata`Q;
    O := Ldata`O;
    basis_L := Ldata`basis_L;
    disc_grp := Ldata`disc_grp;
    to_disc := Ldata`to_disc;

    scale := ScaleForSchofer(d,D,N);
    if Type(Lambda) eq BoolElt then 
        lambda := ElementOfNorm(Q, -d,  O, basis_L);
    else
        lambda := Lambda;
    end if;

    fs := [qExpansionAtoo(eta,1) : eta in etas];
    fs_0 := [qExpansionAt0(eta,1) : eta in etas];

    // Taking care of the principal part at infinity
    log_coeffs := SchoferFormula(fs, d, Q, lambda, scale);

    // Taking care of the principal part at zero
    M := IsOdd(D*N) select 4*D*N else 2*D*N;
    log_coeffs_0 := SchoferFormula0(fs_0, d, Q, lambda, scale, M, disc_grp, to_disc);

    // summing up
    for i->s in log_coeffs do
        log_coeffs[i] +:= log_coeffs_0[i];
    end for;

    return log_coeffs;
end intrinsic;

intrinsic SchoferFormula(eta::EtaQuot, d::RngIntElt, D::RngIntElt, N::RngIntElt, Ldata::QuaternionLatticeData : Lambda := false) -> LogSm
{Return the log of the absolute value of Psi_F_f at the CM point with CM d.}
    return SchoferFormula([eta], d, D, N, Ldata : Lambda := Lambda)[1];
end intrinsic;

intrinsic AbsoluteValuesAtRationalCMPoint(fs::SeqEnum[EtaQuot], d::RngIntElt, Xstar::ShimuraQuot, Ldata::QuaternionLatticeData : Lambda := false) -> SeqEnum[LogSm]
{Returns the absolute value of f for every f in fs at the rational CM point with CM d.}
    vals := [LogSum() : f in fs];
    for i->f in fs do
        div_f := DivisorOfBorcherdsForm(f, Xstar);
        in_support := exists(pt){pt : pt in div_f | pt[1] eq d};
        if in_support then
            if pt[2] lt 0 then vals[i] := LogSum(Infinity()); end if;
            if pt[2] gt 0 then vals[i] := LogSum(0); end if;
        end if;
    end for;
    rest_idxs := [i : i in [1..#fs] | vals[i] eq LogSum()];
    if IsEmpty(rest_idxs) then return vals; end if;
    rest_fs := [fs[i] : i in rest_idxs];
    log_coeffs := SchoferFormula(rest_fs, d, Xstar`D, Xstar`N, Ldata : Lambda := Lambda);
    for i->log_coeff in log_coeffs do
        vals[rest_idxs[i]] := log_coeff;
    end for;
    return vals;
end intrinsic;

intrinsic CandidateDiscriminants(Xstar::ShimuraQuot, curves::SeqEnum[ShimuraQuot] : Exclude := {}, bd := 4) -> SeqEnum
{Returns list of candidate discriminats for Schofer's formula} //'
    rat_pts, quad_pts := RationalandQuadraticCMPoints(Xstar : Exclude := Exclude, coprime_to_level := true, bd := bd);
    return [rat_pts, quad_pts];
end intrinsic;


intrinsic AbsoluteValuesAtCMPoints(Xstar::ShimuraQuot, curves::SeqEnum[ShimuraQuot], all_cm_pts ::SeqEnum, fs::Assoc : MaxNum := 7, Prec := 100,  Exclude := {}, Include := {}) -> SchoferTable, SeqEnum
{Returns the absolute values of y^2 for all degree 2 covers and two hauptmodules at CM points.}
    
    cm_pts_rat := [a : a in all_cm_pts[1]| a[1] notin Exclude];
    cm_pts_quad := [a : a in  all_cm_pts[2] | a[1] notin Exclude];

    keys_fs := [k : k in Keys(fs)];
    all_fs := [fs[k] : k in keys_fs];
   
    if #Include gt 0 then
        include_bd := Maximum([ ClassNumber(d) : d in Include]);
    else
        include_bd := 0;
    end if;

    cm_pts_must_rational := [p : p in cm_pts_rat | p[1] in Include];
    cm_pts_must_quad := [p : p in cm_pts_quad | p[1] in Include];
    other_cm_rat := [p : p in cm_pts_rat | p[1] notin Include];
    other_cm_quad := [p : p in cm_pts_quad | p[1] notin Include];
    need := MaxNum - #Include;
    pt_list_quad := [];
    if #other_cm_rat ge need then
    //need to make space for include points, but otherwise fill up with rational points as much as possible
        pt_list_rat := cm_pts_must_rational cat other_cm_rat[1..need];
        pt_list_quad cat:= cm_pts_must_quad;
    else
        need := need - #other_cm_rat; //update how many we need, after adding rational points
        vprintf ShimuraQuotients, 3: "Still need %o rational points\n", need;
        pt_list_rat := cm_pts_must_rational cat other_cm_rat; //now go search for more points
        Exclude := Exclude join {pt[1] : pt in pt_list_rat};
        bd := Maximum(include_bd*2, 8); //go up to 8 from 4
        new_rat_cm, new_quad_cm := RationalandQuadraticCMPoints(Xstar : bd := bd, Exclude := Exclude, coprime_to_level := true);
        pt_list_rat := pt_list_rat cat new_rat_cm;
        need := need - #new_rat_cm;
        if need gt 0 then
        //now add quadratic points
            other_cm_quad := other_cm_quad cat new_quad_cm;
            if #other_cm_quad ge need then
                pt_list_quad := pt_list_quad cat other_cm_quad[1..need];
            else
                error "Could not find enough points, sorry!";
            end if;
        end if;
    end if;

    table := [[] : f in all_fs];
    // _,_,_,_,_,Q,O,basis_L := ShimuraCurveLattice(Xstar`D,Xstar`N);
    Ldata := ShimuraCurveLattice(Xstar`D,Xstar`N);
    Q := Ldata`Q;
    O := Ldata`O;
    basis_L := Ldata`basis_L;

    vprintf ShimuraQuotients, 2: "\tfinding elements of norm for %o CM point(s)...", #pt_list_rat + #pt_list_quad;
    tt := Realtime();
    lambdas := ElementsOfNorm(Q, [-pt[1] : pt in pt_list_rat cat pt_list_quad], O, basis_L);
    vprintf ShimuraQuotients, 2: " done (%os).\n", Realtime() - tt;
    for j->pt in pt_list_rat do
        d := pt[1];
        tt := Realtime();
        vprintf ShimuraQuotients, 2: "\tabsolute values at rational CM point %o/%o (d = %o)...", j, #pt_list_rat, d;
        vals := AbsoluteValuesAtRationalCMPoint(all_fs, d, Xstar, Ldata : Lambda := lambdas[-d]);
        for i->v in vals do
            Append(~table[i], vals[i]);
        end for;
        vprintf ShimuraQuotients, 2: " done (%os).\n", Realtime() - tt;
    end for;

    for j->pt in pt_list_quad do
        d := pt[1];
        tt := Realtime();
        vprintf ShimuraQuotients, 2: "\tabsolute values at quadratic CM point %o/%o (d = %o)...", j, #pt_list_quad, d;
        norm_val := AbsoluteValuesAtRationalCMPoint(all_fs, d, Xstar, Ldata : Lambda := lambdas[-d]);
        for i->v in norm_val do
            Append(~table[i], norm_val[i]);
        end for;
        vprintf ShimuraQuotients, 2: " done (%os).\n", Realtime() - tt;
    end for;

    all_cm_pts := [cm_pts_rat, cm_pts_quad];
    ds := [pt[1] : pt in pt_list_rat] cat [ pt[1] : pt in pt_list_quad ];

    schofer_tab := CreateSchoferTable(table, keys_fs, ds, curves,Xstar);
    schofer_tab`BorcherdsForms := fs;

    return schofer_tab, all_cm_pts;
end intrinsic;

function find_degs(abs_schofer_tab)
//returns list of degrees of fields of definition for each d, on the star curve
    ds := abs_schofer_tab`Discs;
    fldsofdef := abs_schofer_tab`FldsOfDefn;
    cid := abs_schofer_tab`Xstar`CurveID;
    degs := [];
    for d in ds do
        Append(~degs, Degree(fldsofdef[cid][d][1]));
    end for;
    return degs;
end function;

function find_signs_hauptmodul(s, stilde, ds, degs)
    inf_zero_indices := [Index(s,0), Index(stilde,0), Index(s,Infinity())];
    assert stilde[inf_zero_indices[3]] eq Infinity();
    scale_tilde := stilde[Index(s,0)];
    scale := s[Index(stilde,0)];

    rat_idxs := [i : i in [1..#s] | i notin inf_zero_indices and degs[i] eq 1];
    signs := &cat[[[eps1, eps2] : eps1,eps2 in [-1,1] | eps1*s[i]/scale + eps2*stilde[i]/scale_tilde eq 1] : i in rat_idxs];
    s_new := [* ss/scale^degs[i] : i->ss in s *];
    stilde_new := [* sstilde/scale_tilde^degs[i] : i->sstilde in stilde *];
    for j->idx in rat_idxs do
        s_new[idx] := signs[j][1]*s_new[idx];
        stilde_new[idx] := signs[j][2]*stilde_new[idx];
    end for;
    return s_new, stilde_new, scale, scale_tilde;
end function;


// This is following [GR, Section 5]
intrinsic FieldsOfDefinitionOfCMPoint(X::ShimuraQuot, d::RngIntElt) -> List
{Return possible fields of definition for CM point with CM by d on X.}
    // require IsFundamentalDiscriminant(d) : "Field of definition currently only supports maximal orders";
    R := QuadraticOrder(BinaryQuadraticForms(d));
    K := NumberField(R);
    f := Conductor(R);
    H_R := RingClassField(R);
    D := X`D;
    N := X`N;
    D_R := &*[Integers()| p : p in PrimeDivisors(D) | KroneckerCharacter(d)(p) eq -1];
    N_R := &*[Integers()| p : p in PrimeDivisors(N) | KroneckerCharacter(d)(p) eq 1 or (f mod p eq 0)];   
    N_star_R := &*[Integers()| p : p in PrimeDivisors(N) | (KroneckerCharacter(d)(p) eq 1) and (f mod p ne 0)];
    assert GCD(D_R * N_star_R, Discriminant(R)) eq 1;
    assert GCD(D_R*N_R, Discriminant(R)) eq GCD(N,f);

    // Proposition 5.6
    if (Discriminant(R) mod ((D*N) div (D_R*N_star_R))) ne 0 then
        return [* *];
    end if;

    rec := ArtinMap(H_R);

    // also number of points is 2^PrimeDivisors(D_R*N_R) * ClassNumber(R)

    // Theorem 5.8 - Shimura reciprocity
    
    // setting up number fields
    H_R_NF := NumberField(H_R);
    abs_H_R := AbsoluteField(H_R_NF);
    _, H_R_to_abs := IsIsomorphic(H_R_NF, abs_H_R);

    // setting up Picard groups 
    PicR, mPicR := PicardGroup(R);
    A, PicR_to_A := PicR / (2*PicR);
    B := QuaternionAlgebra(D);

    al_action := AssociativeArray();

    // Theorem 5.12 (1) and Lemma 5.10 for complex conjugation
    m := D_R*N_star_R;
    
    frakas := [];
    for a in A do
        alift := a@@PicR_to_A;
        // circumventing a bug in Magma in mPicR
        if (alift eq PicR!0) then
            fraka := 1*R;
        else
            fraka := mPicR(alift);
        end if;
        B_fraka := QuaternionAlgebra(Rationals(), d, m*Norm(fraka));
        if IsIsomorphic(B_fraka, B) then
            Append(~frakas, fraka);
            if IsFundamentalDiscriminant(d) then break; end if; // see [GR], Remark 5.11
        end if;
    end for;
    require #frakas gt 0 : "Error in field of definition - could not find a fractional ideal for complex conjugation!";
    sigma_as := [rec(fraka) : fraka in frakas];
    abs_sig_as := [H_R_to_abs^(-1)*sigma_a*H_R_to_abs : sigma_a in sigma_as];

    // Lemma 5.9
    fixed_sub_gens := [];
    unknown_quotients := 0;
    als_DN := [Q : Q in Divisors(D*N) | GCD(Q, (D*N) div Q) eq 1];
    // we already have m so use mm just to be safe
    for mm in als_DN do
        al_is_gal := ((D*N) div (D_R*N_R)) mod mm eq 0;
        if al_is_gal then
            frakb := &*[Parent(1*Integers(K)) | pa[1]^(pa[2] div 2) : pa in Factorization(mm*Integers(K))];
            assert Norm(frakb) eq mm;
            al_action[mm] := H_R_to_abs^(-1)*rec(frakb)*H_R_to_abs;
        end if;
        // !! TODO : figure out what to do if it is not Galois
    end for;

    known_al := Keys(al_action);
    S := Subsets(known_al);
    for s in S do
        if #s eq 0 then
            al_action[1] := hom<abs_H_R->abs_H_R | abs_H_R.1>;
        else
            prod := 1;
            for w in s do
                prev_prod := prod;
                prod := AtkinLehnerMul(w,prod,D*N);
                if not IsDefined(al_action, prod) then
                    al_action[prod] := al_action[prev_prod]*al_action[w];
                end if;
            end for;
        end if;
    end for;

    fixed_by := [al_action[mm] : mm in X`W meet Keys(al_action)];

    Q_P := FixedField(abs_H_R, fixed_by);
    Q_Ps := [* Q_P *];

    // Handle complex conjugation 
    has_cc, cc := HasComplexConjugate(abs_H_R);
    if not has_cc then
        gal, auts, gal_to_auts := AutomorphismGroup(abs_H_R);
        // elements that restrict to the complex conjugation on K
        cc_candidates := [g : g in gal | Order(g) eq 2 and gal_to_auts(g)(K.1) eq ComplexConjugate(K.1)];  
        cc_reps := [gal_to_auts(cc) : cc in cc_candidates];
    else
        cc_reps := [cc];
    end if;
    sigmas := [hom<abs_H_R -> abs_H_R | cc(abs_sig_a(abs_H_R.1))> : abs_sig_a in abs_sig_as, cc in cc_reps];
    if m eq 1 then 
        return [* FixedField(Q_P, [sigma]) : Q_P in Q_Ps, sigma in sigmas *];
    end if;
    Q_Ps := [* *];
    known_al := Keys(al_action);
    for sigma in sigmas do
        al_action[m] := sigma;
        for w in known_al do
            prod := AtkinLehnerMul(w,m,D*N);
            al_action[prod] := al_action[m]*al_action[w];
        end for;
        fixed_by := [al_action[mm] : mm in X`W meet Keys(al_action)];
        Q_P := FixedField(abs_H_R, fixed_by);
        Append(~Q_Ps, Q_P);
    end for;

    return Q_Ps;

end intrinsic;

intrinsic FieldsOfDefinitionOfCMPointFast(X::ShimuraQuot, d::RngIntElt) -> List
{Faster variant of FieldsOfDefinitionOfCMPoint: returns the possible fields of
 definition of the CM point with CM by d on X, built via Magma's AbelianExtension
 inside the (smaller) Atkin-Lehner-fixed field A_abs rather than the full ring class
 field H_R + ArtinMap(H_R).  Returns the same set of fields (up to isomorphism) as
 FieldsOfDefinitionOfCMPoint.  See arXiv:math/0612732v2, Appendix.}
    D := X`D;
    N := X`N;
    W := X`W;
    R := QuadraticOrder(BinaryQuadraticForms(d));
    K := NumberField(R);
    f := Conductor(R);
    OK := Integers(K);
    chi := KroneckerCharacter(d);

    D_R      := &*[Integers()| p : p in PrimeDivisors(D) | chi(p) eq -1];
    N_R      := &*[Integers()| p : p in PrimeDivisors(N) | chi(p) eq 1 or (f mod p eq 0)];
    N_star_R := &*[Integers()| p : p in PrimeDivisors(N) | chi(p) eq 1 and (f mod p ne 0)];

    // Proposition 5.6
    if (Discriminant(R) mod ((D*N) div (D_R*N_star_R))) ne 0 then
        return [* *];
    end if;

    mG := NormGroup(RingClassField(R));   // Gal(H_R/K) -> ideals of O_K
    G  := Domain(mG);

    // Atkin-Lehner-Galois subgroup alSub_W <= G, restricted to the quotient W
    // (Lemma 5.9): generated by the Artin classes of the "square root" ideals frakb_mm
    // for the al_is_gal operators w_mm with mm in W.
    al_gen := [];
    for mm in W do
        if mm eq 1 then continue; end if;
        if ((D*N) div (D_R*N_R)) mod mm ne 0 then continue; end if;   // al_is_gal
        frakb := &*[Parent(1*OK) | pa[1]^(pa[2] div 2) : pa in Factorization(mm*OK)];
        Append(~al_gen, frakb @@ mG);
    end for;

    // A_abs = H_R fixed by alSub_W, as an abelian extension of K.
    _, pi_quo := quo<G | sub<G | al_gen>>;
    Aext := AbelianExtension(Inverse(pi_quo)*mG);
    NFA  := NumberField(Aext);
    Aabs := AbsoluteField(NFA);
    _, NFA_to_abs := IsIsomorphic(NFA, Aabs);
    KinA := Aabs!K.1;

    // Is complex conjugation active on this quotient?  Lemma (CC): complex conjugation
    // is realised by w_m . sigma_a with m = D_R*N_star_R; it identifies P with bar P on
    // the quotient iff some mm in W is a reflection m * w0 with w0 al_is_gal.  In that
    // case the reflection restricts to A_abs as c . sigma_a . sigma_{w0}, so we must
    // include the extra Atkin-Lehner twist sigma_{w0} (trivial only when w0 = 1, e.g.
    // for the full quotient).
    m := D_R*N_star_R;
    cc_active := false;
    w0 := 1;
    for mm in W do
        w := AtkinLehnerMul(m, mm, D*N);
        if ((D*N) div (D_R*N_R)) mod w eq 0 then   // w is al_is_gal
            cc_active := true;
            w0 := w;
            break;
        end if;
    end for;

    if not cc_active then
        return [* Aabs *];   // field of definition contains K
    end if;

    // Complex-conjugation candidate(s).  The slow function uses HasComplexConjugate on the
    // FULL ring class field H_R: a single canonical conjugation when H_R is CM, else every
    // order-2 automorphism restricting to complex conjugation on K.  H_R is CM iff complex
    // conjugation is central iff cc acts trivially by inversion on Gal(H_R/K) = Pic(R), i.e.
    // iff Pic(R) has exponent <= 2 (here G = Domain(mG) ~ Pic(R)).  We reproduce that on the
    // smaller A_abs (whose cc-candidates restrict from those of H_R, giving the same fixed
    // fields):
    //   exponent <= 2  ->  H_R CM  ->  the unique canonical conjugation (one field per a);
    //   exponent  > 2  ->  enumerate all cc-candidates (multiple genuine fields).
    if Exponent(G) le 2 then
        has_cc, cc := HasComplexConjugate(Aabs);   // A_abs (subfield of CM H_R, contains K) is CM
        assert has_cc;
    else
        // A_abs is not CM (Pic(R) has exponent > 2), so HasComplexConjugate fails -- but complex
        // conjugation is still ONE automorphism (A_abs/Q is Galois): the unique one sending A_abs.1
        // to the complex-conjugate root of its defining polynomial.  GR (Lemma CC / Thm mainCM) use
        // this c abstractly and pin only sigma_a; we build c directly here rather than enumerating
        // every order-2 reflection via AutomorphismGroup (slow, and over-generating -- it returns
        // all reflections c.sigma, of which only this one is genuinely complex conjugation).
        f := DefiningPolynomial(Aabs);
        prec := 40 + 4*Degree(Aabs);
        target := ComplexConjugate(Conjugates(Aabs.1 : Precision := prec)[1]);
        roots_A := [r[1] : r in Roots(f, Aabs)];
        diffs := [Abs(Conjugates(r : Precision := prec)[1] - target) : r in roots_A];
        _, idx := Minimum(diffs);
        cc := hom<Aabs -> Aabs | roots_A[idx]>;
    end if;
    ccs := [* cc *];

    // Valid classes [a] in Pic(R)/Pic(R)^2 with B_D ~ (-s, m*N(a))_Q  (Lemma CC /
    // Remark uniqueness); their Artin automorphisms sigma_a restricted to A_abs.
    B := QuaternionAlgebra(D);
    PicR, mPicR := PicardGroup(R);
    A2, PicR_to_A2 := PicR / (2*PicR);
    am := ArtinMap(Aext);

    // Extra Atkin-Lehner twist sigma_{w0} on A_abs from the reflection m*w0 in W.
    frakb_w0 := &*[Parent(1*OK) | pa[1]^(pa[2] div 2) : pa in Factorization(w0*OK)];
    tw := (NFA_to_abs^-1)*am(frakb_w0)*NFA_to_abs;

    sigma_as := [* *];
    for a in A2 do
        alift := a@@PicR_to_A2;
        fraka := (alift eq PicR!0) select 1*R else mPicR(alift);
        if not IsIsomorphic(QuaternionAlgebra(Rationals(), d, m*Norm(fraka)), B) then
            continue;
        end if;
        Append(~sigma_as, (NFA_to_abs^-1)*am(OK !! fraka)*NFA_to_abs);
        if IsFundamentalDiscriminant(d) then break; end if;   // [GR] Remark 5.11
    end for;
    require #sigma_as gt 0 :
        "Error in field of definition - could not find a fractional ideal for complex conjugation!";

    // Q(P) = fixed field of the reflection c . sigma_a . sigma_{w0}; collect distinct.
    fields := [* *];
    for c in ccs do
        for sa in sigma_as do
            sigma := hom<Aabs -> Aabs | tw(c(sa(Aabs.1)))>;
            F := FixedField(Aabs, [sigma]);
            is_new := true;
            for FF in fields do
                if IsIsomorphic(F, FF) then is_new := false; break; end if;
            end for;
            if is_new then Append(~fields, F); end if;
        end for;
    end for;

    return fields;

end intrinsic;

procedure replace_column(schofer_tab, d, dnew, is_log)
    //add column associated to dnew remove column associated to d
    table := schofer_tab`Values;
    ds := schofer_tab`Discs;
    Xstar := schofer_tab`Xstar;
    curveid := Xstar`CurveID;
    fs := schofer_tab`BorcherdsForms;
    keys_fs := schofer_tab`Keys_fs;
    row_scales := schofer_tab`RowScales;
    all_fs := [fs[k] : k in keys_fs];
    UpdateFieldsOfDefn(schofer_tab, dnew);
    flds := (schofer_tab`FldsOfDefn)[curveid][dnew];
    assert #flds eq 1;
    deg := Degree(flds[1]);
    d_idx := Index(ds,d);
    ds[d_idx] := dnew;
    Ldata := ShimuraCurveLattice(Xstar`D,Xstar`N);
    norm_val := AbsoluteValuesAtRationalCMPoint(all_fs, dnew, Xstar, Ldata);
    for i->v in norm_val do
        // table[i][d_idx] := norm_val[i]/row_scales[i]^deg;
        if is_log then
            table[i][d_idx] := norm_val[i]-deg*row_scales[i];
        else
            table[i][d_idx] := RationalNumber(norm_val[i]-deg*row_scales[i]);
        end if;
    end for;
    schofer_tab`Values := table;
    schofer_tab`Discs := ds;
    curves := schofer_tab`Curves;
    return;
end procedure;


function find_y2_scales(schofer_table)
    ds := schofer_table`Discs;
    degs := find_degs(schofer_table);
    ratds := [d : i->d in ds | degs[i] eq 1];
    table := schofer_table`Values;
    keys_fs := schofer_table`Keys_fs;
    k_idxs := schofer_table`K_idxs;
    fldsofdef := schofer_table`FldsOfDefn;

    //Scale the y2 rows of the table

    scale_factors :=[];
    for i in k_idxs do
        if exists(j1){j : j->d1 in ratds  | #fldsofdef[keys_fs[i]][d1] eq 1 and Degree(fldsofdef[keys_fs[i]][d1][1]) eq 1 and table[i][j] ne LogSum(Infinity()) and table[i][j] ne LogSum(0)} then
            //then we have a rational point on X
            // d1 := ratds[j1];
            v1 := table[i][j1];
            // scale, _ := SquareFree(v1);
            log_scale := SquareFree(v1);
            // Append(~scale_factors, AbsoluteValue(scale));
            Append(~scale_factors, log_scale);
        else 
            found_j1 := exists(j1){j : j->d1 in ratds  | #fldsofdef[keys_fs[i]][d1] le 2 and {Degree(fldsofdef[keys_fs[i]][d1][k]) : k in [1..#fldsofdef[keys_fs[i]][d1]]} subset {1,2} and table[i][j] ne LogSum(Infinity()) and table[i][j] ne LogSum(0)};
            found_j2 := found_j1 and exists(j2){j : j->d2 in ratds  | #fldsofdef[keys_fs[i]][d2] le 2 and {Degree(fldsofdef[keys_fs[i]][d2][k]) : k in [1..#fldsofdef[keys_fs[i]][d2]]} subset {1,2}  and table[i][j] ne LogSum(Infinity()) and ratds[j1] ne d2 and table[i][j] ne LogSum(0)};
            // Graceful: without two suitable rational CM points we cannot pin this cover's y2-scale.
            // Leave the row unscaled (placeholder); the cover's constraints then come out inconsistent
            // in EquationsOfCovers, so it is deferred and (if a parent is computed) recovered as a quotient.
            if not found_j2 then
                Append(~scale_factors, LogSum(1));
                vprintf ShimuraQuotients, 1 : "  Could not pin y2-scale (sparse CM data); leaving a cover unscaled to be deferred downstream.\n";
                continue;
            end if;
            //otherwise we find two points that are potentially over quadratic fields
            v1 := table[i][j1];
            v2 := table[i][j2];
            d1 := ratds[j1];
            d2 := ratds[j2];
            quad_idx_1 := 1;
            quad_idx_2 := 1;
            if fldsofdef[keys_fs[i]][d1][1] eq Rationals() then
                quad_idx_1 := 2;
            end if;
            if fldsofdef[keys_fs[i]][d2][1] eq Rationals() then
                quad_idx_2 := 2;
            end if;
            d1 := Discriminant(MaximalOrder(fldsofdef[keys_fs[i]][d1][quad_idx_1]));
            d2 := Discriminant(MaximalOrder(fldsofdef[keys_fs[i]][d2][quad_idx_2]));
            log_scale1 := SquareFree(v1 - LogSum(AbsoluteValue(d1)));
            log_scale2 := SquareFree(v1);
            // scale1, _ := SquareFree(v1/d1); //two possibilities
            // scale2, _ := SquareFree(v1);
            // if IsSquare(scale1*v2/d2) then
            if IsSquare( v2 - LogSum(AbsoluteValue(d2)) - log_scale1) then
                // Append(~scale_factors, AbsoluteValue(scale1));
                Append(~scale_factors, log_scale1);
            elif IsSquare(v2 - log_scale2) then
                // Append(~scale_factors, AbsoluteValue(scale2));
                Append(~scale_factors, log_scale2);
            else
                Append(~scale_factors, LogSum(1));
                vprintf ShimuraQuotients, 1 : "  y2-scale IsSquare check failed for a cover; leaving it unscaled to be deferred downstream.\n";
            end if;
        end if;
    end for;
    return scale_factors;

end function;

function find_y2_signs(abs_schofer_tab)
    //find signs of y^2 for rational CM point d on each y^2
    //in keys_fs, where d_idx is index of column of d in table
    table := abs_schofer_tab`Values;
    keys_fs := abs_schofer_tab`Keys_fs;
    k_idxs := abs_schofer_tab`K_idxs;
    flds := abs_schofer_tab`FldsOfDefn;
    ds := abs_schofer_tab`Discs;
    degs := find_degs(abs_schofer_tab);
    rat_idxs := [i : i in [1..#ds] | degs[i] eq 1];
    quad_idxs := [i : i in [1..#ds] | degs[i] eq 2];

    for d_idx->d in ds do
        for k->i in k_idxs do
            if table[i][d_idx] eq Infinity() then continue; end if;
            if table[i][d_idx] eq 0 then continue; end if;
            fields := flds[keys_fs[i]][d];
            possible_answers := [* *];
            for eps in [-1,1] do
                y2 := eps*table[i][d_idx];
                for F in fields do
                    is_sqr, y := IsSquare(F!y2);
                    if (is_sqr) then
                        if d_idx in rat_idxs then
                            if (Type(F) eq FldRat) or (Degree(F) eq Degree(sub<F|y>)) then
                                Append(~possible_answers, <F,eps,y>);
                            end if;
                        elif d_idx in quad_idxs then
                            if Degree(MinimalPolynomial(y)) eq 1 then
                                if NormEquation(F, y2) then
                                    Append(~possible_answers, <F,eps,y>); //this is just a norm value
                                end if;
                            end if;
                        end if;
                    end if;
                end for;
            end for;
            // When the sign is determined (cover value of degree <= 2, i.e. f_q(s) a square in the
            // star field), apply it.  Otherwise -- the cover value is degree 4 (f_q(s) NOT a square,
            // which happens for the exponent>2 quadratic discs) -- the sign is neither determinable
            // from the y^2 norm alone nor needed: the cover equations f_q are fixed by the rational
            // and degree-2 cover points, and the special-fiber method uses only the star Hauptmodul
            // s-value of such a disc.  Leave |y^2| untouched and skip.
            if #possible_answers eq 1 then
                eps := possible_answers[1][2];
                table[i][d_idx] :=  eps * table[i][d_idx];
            end if;
        end for;
    end for;
    return table;
end function;

function find_hauptmodul_signs_quadratic(abs_schofer_tab, d, d_idx)
    //find signs of hauptmodul for quadratic CM point d on each hauptmodul
    //in keys_fs, where d_idx is index of column of d in table
    table := abs_schofer_tab`Values;
    keys_fs := abs_schofer_tab`Keys_fs;
    k_idxs := abs_schofer_tab`K_idxs;
    s_idx := abs_schofer_tab`sIndex;
    stilde_idx := abs_schofer_tab`sTildeIndex;
    flds := abs_schofer_tab`FldsOfDefn;

    for k->i in k_idxs do
        if table[i][d_idx] eq Infinity() then continue; end if;
        if table[i][d_idx] eq 0 then continue; end if;
        K := flds[keys_fs[i]][d][1];
        norm_s := table[s_idx][d_idx];
        norm_stilde := table[stilde_idx][d_idx];
         _<x> := PolynomialRing(Rationals());
        signs := [[1,1], [1,-1],[-1,1],[-1,-1]];
        minpolys := [];
        for eps in signs do
                trace := 1 - eps[1]*norm_stilde +  eps[2]*norm_s;
                Append(~minpolys, x^2 - trace*x + eps[2]*norm_s);
        end for;
        roots := [Roots(p,K) : p in minpolys];
        good_inds := [i : i->r in roots | #r ne 0 and not(&and[rt[1] in Rationals() : rt in r])];
        if #good_inds eq 1 then
            table[s_idx][d_idx] := minpolys[good_inds[1]];
            norm_s := Coefficient(minpolys[good_inds[1]], 0);
            trace_s := - Coefficient(minpolys[good_inds[1]], 1);
            table[stilde_idx][d_idx] := x^2 - (2- trace_s)*x + (1- trace_s + norm_s);
            return true, table;
        else
            vprintf ShimuraQuotients, 3: "We need that there is a unique minpoly left after filtering by roots, but we found %o good indices\n", #good_inds;
            if #good_inds eq 0 then
                error "No good indices found after filtering by roots";
            end if;
            return false, _;
        end if;
    end for;
end function;

intrinsic ValuesAtCMPoints(abs_schofer_tab::SchoferTable, all_cm_pts::SeqEnum : Exclude := {}) -> SchoferTable
    {}
    allds := abs_schofer_tab`Discs;
    table := abs_schofer_tab`Values;
    Xstar := abs_schofer_tab`Xstar;
    row_scales := abs_schofer_tab`RowScales;

    //scale the hauptmodul rows
    s_idx := abs_schofer_tab`sIndex;
    stilde_idx := abs_schofer_tab`sTildeIndex;
    s := table[s_idx];
    stilde := table[stilde_idx];

    table :=[* [* x : i->x in t *] : t in table *];
    scale_tilde := stilde[Index(s,LogSum(0))];
    scale := s[Index(stilde,LogSum(0))];
   
    row_scales[s_idx] +:= scale;
    row_scales[stilde_idx] +:= scale_tilde;
   
    k_idxs := abs_schofer_tab`K_idxs;
    abs_schofer_tab`Values := table;

    //Scale the y2 rows of the table
    scale_factors := find_y2_scales(abs_schofer_tab);

    degs := find_degs(abs_schofer_tab);
    for i->k in k_idxs do
        for j in [1 .. #table[k]] do
            // table[k][j] := table[k][j]/scale_factors[i]^degs[j];
            table[k][j] -:= degs[j]*scale_factors[i];
        end for;
        // row_scales[k] := row_scales[k]*scale_factors[i];
        row_scales[k] +:= scale_factors[i];
    end for;

    // make table values into rational numbers
    abs_schofer_tab`Values := [*[*RationalNumber(x) : x in y*] : y in table*];
    abs_schofer_tab`RowScales := row_scales;
    
    //find signs on the hauptmodul rows
    table := abs_schofer_tab`Values;
    s := table[s_idx];
    stilde := table[stilde_idx];
    degs := find_degs(abs_schofer_tab);
    s, stilde := find_signs_hauptmodul(s, stilde, allds, degs);
    table[s_idx] := s;
    table[stilde_idx] := stilde;
    abs_schofer_tab`Values := table;

    quad_idxs := [i : i in [1..#allds] | degs[i] eq 2];
    used_ds := Set(allds);
    for i in [1..#allds] do
        if i notin quad_idxs then continue; end if; //only do quadratic points
        currd := allds[i];
        success, new_table := find_hauptmodul_signs_quadratic(abs_schofer_tab, currd, i);
        while not success do
            vprintf ShimuraQuotients, 1: "We need that there is a unique minpoly left after filtering by roots so we are replacing %o.\n", currd;
            Include(~used_ds, currd);
            candidates := Set([pt[1] : pt in all_cm_pts[2]]) diff used_ds diff Exclude;
            if #candidates eq 0 then
                error "No possible choices of CM points left which we can pin down the correct minpoly";
            end if;
            newd := Reverse(Sort(SetToSequence(candidates)))[1];
            vprintf ShimuraQuotients, 1: "Replacing %o with %o\n", currd, newd;
            replace_column(abs_schofer_tab, currd, newd, false);
            Include(~used_ds, newd);
            currd := newd;
            success, new_table := find_hauptmodul_signs_quadratic(abs_schofer_tab, currd, i);
        end while;
        abs_schofer_tab`Values := new_table;
    end for;

    table := find_y2_signs(abs_schofer_tab);

    schofer_table := CreateSchoferTable(table, abs_schofer_tab`Keys_fs, abs_schofer_tab`Discs, abs_schofer_tab`Curves, Xstar);
    return schofer_table;
end intrinsic;

intrinsic ReduceTable(schofer_tab::SchoferTable)
    {}
    table := schofer_tab`Values;
    allds := schofer_tab`Discs;
    degs := find_degs(schofer_tab);
    scales := [];
    rat_idxs := [i : i in [1..#degs] | degs[i] eq 1];
    for t in table do
        // xs := [x : i->x in t | i le num_rat_ds and x notin [0, Infinity()] ];
        xs := [x : i->x in t | i in rat_idxs and x notin [LogSum(0), LogSum(Infinity())] ];
        //xs := [x : x in t | x notin [0, Infinity()]];
        /*
        ps := &join[Set(PrimeDivisors(Numerator(x))) : x in xs];
        ps join:= &join[Set(PrimeDivisors(Denominator(x))) : x in xs];
        ps := [p : p in ps];
        */
        ps := [p : p in &join[Keys(x`log_coeffs) : x in xs]];
        // vals := [[Valuation(x,p) : x in xs ] : p in ps];
        vals := [[IsDefined(x`log_coeffs, p) select x`log_coeffs[p] else 0 : x in xs ] : p in ps];
        mins := [Minimum([<AbsoluteValue(v),v> : v in valp]) : valp in vals];
        scale := &+([LogSum()] cat [LogSum(mins[i][2], p) : i->p in ps]);
        Append(~scales, scale);
    end for;
    schofer_tab`Values :=  [[x - degs[j]*scales[i] : j->x in t] : i->t in table ];
    schofer_tab`RowScales := scales;
    return;
end intrinsic;

intrinsic ValuesAtCMPoints(Xstar::ShimuraQuot, curves::SeqEnum[ShimuraQuot] : MaxNum := 7, Prec := 100, Exclude := {}, Include := {}) -> SeqEnum, SeqEnum, SeqEnum
{Returns the values of y^2 for all degree 2 covers and two hauptmodules at CM points.}
    fs := BorcherdsForms(Xstar, curves : Prec := Prec);
    d_divs := &cat[[T[1]: T in  DivisorOfBorcherdsForm(f, Xstar)] : f in [fs[-1], fs[-2]]]; //include zero infinity of hauptmoduls
    all_cm_pts := CandidateDiscriminants(Xstar, curves);
    abs_schofer_tab, all_cm_pts := AbsoluteValuesAtCMPoints(Xstar, curves, all_cm_pts, fs : MaxNum := MaxNum, Prec := Prec, Exclude := {}, Include := Set(d_divs));
    ReduceTable(abs_schofer_tab);
    schofer_tab := ValuesAtCMPoints(abs_schofer_tab, all_cm_pts : Exclude := Exclude);
    return schofer_tab;
end intrinsic;

/* Bibliography
[GR] - Gonzales, Rotger - non-elliptic shimura curves of genus one
*/
