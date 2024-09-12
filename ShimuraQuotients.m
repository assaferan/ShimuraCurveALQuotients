// Lower and upper bounds for number of points in reduction mod p, 
// scaled by 12/p-1, see [HH]

function Omega(n)
    return #PrimeDivisors(n);
end function;

function LowerBound(D, N)
    ps := PrimeDivisors(N);
    c_D := EulerPhi(D);
    c_N := N * &*[Rationals() | 1 + 1/p : p in ps];
    return c_D * c_N / (2^Omega(D*N));
end function;

function UpperBound(p)
    return 24*(1+p^2)/(p-1);
end function;

function MinimalNumberOfALinQuotient(D, N)
    ps := PrimeDivisors(N);
    c_D := EulerPhi(D);
    c_N := N * &*[Rationals() | 1 + 1/p : p in ps];
    p := 2;
    while (N mod p eq 0) do
	p := NextPrime(p);
    end while;
    return Maximum(0,Ceiling(Log(2, (c_D * c_N) / UpperBound(p))));
end function;

procedure VerifyBound(: r := 7)
    ps := PrimesUpTo(100);
    assert LowerBound(1, &*ps[1..r]) gt UpperBound(ps[r+1]);
    assert LowerBound(&*ps[1..r], 1) gt UpperBound(ps[r+1]);
end procedure;

function FindMaximalN(r)
    ps := PrimesUpTo(100);
    return Floor(UpperBound(ps[r])*2^(r-1));
end function;

// Using https://math.stackexchange.com/questions/301837/is-the-euler-phi-function-bounded-below EEDDIITT -> TODO: Find a reference
function FindMaximalD(r)
    // Using phi(D) >= Sqrt(D/2)
    // Enough that Sqrt(D/2) gt ub
    // D gt 2*ub^2;
//  return 1290;
    
    ps := PrimesUpTo(100);
    ub := Ceiling(UpperBound(ps[r]));
    return 2*ub^2;
    
    ps := PrimesUpTo(ub+1);
    maxD := 1;
    for num_primes in [1..r] do
	A := CartesianPower(ps, num_primes);
	for a in A do
	    D := &*a;
	    if not IsSquarefree(D) then
		continue;
	    end if;
	    if EulerPhi(D) gt ub then
		break;
	    else
		if D gt maxD then
		    maxD := D;
		end if;
	    end if;
	end for;
    end for;
    return maxD;
end function;

function FindPairs(: r := 7)
    pairs := [];
    N0 := FindMaximalN(r);
    D0 := FindMaximalD(r);
    Ds := [D : D in [1..D0] | IsSquarefree(D) and IsEven(Omega(D))];
    for D in Ds do
	print "D = ", D;
	Nmax := Floor(N0 / EulerPhi(D));
	for N in [1..Nmax] do
	    /*
	    if (N mod 1000 eq 0) then
		print "N =", N;
	    end if;
	   */
	    p := 2;
	    while (N mod p eq 0) do
		p := NextPrime(p);
	    end while;
	    if (LowerBound(D, N) le UpperBound(p)) then
		Append(~pairs, <D, N>);
	    end if;
	end for;
    end for;
    return pairs;
end function;

// 2 corresponds to cusps, 3, 4
function NumberOfEllipticPoints(D, N, order)
    if order eq 2 then
	if D eq 1 then
	    return &+[EulerPhi(GCD(d, N div d)) : d in Divisors(N)];
	else
	    return 0;
	end if;
    end if;
    if order eq 4 then
	Q := 4;
    end if;
    if order eq 3 then 
	Q := 9;
    end if;
    if (N mod Q eq 0) then
	return 0;
    end if;
    primesD := PrimeDivisors(D);
    primesN := PrimeDivisors(N);
    e_D := &*[Integers() | 1 - KroneckerSymbol(-order, p) : p in primesD]; 
    e_N := &*[Integers() | 1 + KroneckerSymbol(-order, p) : p in primesN];
    return e_D * e_N;
end function;

function GenusShimuraCurve(D, N)
    phiD := EulerPhi(D);
    primes := PrimeDivisors(N);
    P1N := Floor(N * &*[Rationals() | 1 + 1/p : p in primes]);
    g := 1 + phiD * P1N / 12;
    e := AssociativeArray();
    for h in [2,3,4] do
	e[h] := NumberOfEllipticPoints(D, N, h);
	g -:= e[h] / h;
    end for;
    assert IsIntegral(g);
    return Floor(g);
end function;
    
function PsiOgg(p, n)
    if (n eq 1) then
	return 1;
    end if;
    is_prime_power, q, k := IsPrimePower(n);
    if (is_prime_power) then
	return ((q eq p) select Floor(p^k *(1 + 1/p)) else 1);
    end if;
    fac := Factorization(n);
    return &*[Integers() | PsiOgg(p, pe[1]^pe[2]) : pe in fac];
end function;

function LegendreSymbol(R, p)
    f := Conductor(R);
    if (f mod p eq 0) then
	return 1;
    end if;
    ZF := MaximalOrder(R);
    return KroneckerCharacter(Discriminant(ZF))(p);
end function;

// Here R is the quadratic order and 
// O is a fixed quaternion Eichler order of level F in 
// the quaternion algebra B of discriminant D. 
// based on Theorem 2 in [Ogg]
function NuOgg(p, R, D, F)
    if (D mod p eq 0) then
	return 1 - LegendreSymbol(R, p);
    end if;
    if Valuation(F, p) eq 1 then
	return 1 + LegendreSymbol(R, p);
    end if;
    assert Valuation(F, p) ge 2;
    f := Conductor(R);
    ZF := MaximalOrder(R);
    chi := KroneckerCharacter(Discriminant(ZF));
    k := Valuation(f, p);
    K := Valuation(F, p);
    if (K ge 2*(1 + k)) then
	if (chi(p) eq 1) then
	    return 2*PsiOgg(p, f);
	end if;
	return 0;
    end if;
    if (K eq 1 + 2*k) then
	if (chi(p) eq 1) then
	    return 2*PsiOgg(p, f);
	end if;
	if (chi(p) eq 0) then
	    return p^k;
	end if;
	assert chi(p) eq -1;
	return 0;
    end if;
    if (K eq 2*k) then
	return p^(k-1)*(p+1+chi(p));
    end if;
    if (K le 2*k - 1) then
	if IsEven(K) then
	    return p^(k div 2) + p^(k div 2 - 1);
	else
	    return 2*p^(k div 2);
	end if;
    end if;
    // Should not reach here
    assert false;
end function;

// Quotient by w_m, m divides DN, following [Ogg]

// The number of the fixed points of w_m on X_0(D,N) 
function NumFixedPoints(D, N, m)
    if (m eq 2) then
	orders := [MaximalOrder(QuadraticField(-1)), 
		   MaximalOrder(QuadraticField(-2))];
    elif (m mod 4 eq 3) then
	F := QuadraticField(-m);
	_, sqrt_minus_m := IsSquare(F!(-m));
	O := MaximalOrder(F);
	alpha := (1 + sqrt_minus_m)/2;
	orders := [sub<O | 1, alpha>, sub<O | 1, 2*alpha>];
    else
	F := QuadraticField(-m);
	_, sqrt_minus_m := IsSquare(F!(-m));
	O := MaximalOrder(F);
	orders := [sub<O | 1, sqrt_minus_m>];
    end if;
    e := 0;
    for R in orders do
	h := PicardNumber(R);
	// Using formula (4) in [Ogg]
	prod := &*[Integers() | 
		  NuOgg(p, R, D, N) : p in PrimeDivisors(D*N) | m mod p ne 0]; 
	e +:= h*prod;
    end for;
    if (D eq 1) and (m eq 4) then
	M := N div 4;
	num_fixed_cusps := &+[Integers() | 
			     EulerPhi(GCD(d, M div d)) : d in Divisors(M)];
	e +:= num_fixed_cusps;
    end if;
    return e;
end function;
    
function GenusShimuraCurveQuotientSingleAL(D, N, m)
    e := NumFixedPoints(D, N, m);
    g_big := GenusShimuraCurve(D, N);
    g := (g_big+1)/2 - e/4;
    assert IsIntegral(g);
    return Floor(g);
end function;

function GenusShimuraCurveQuotient(D, N, als)
    total_e := 0;
    for al in als do
	assert GCD(al, (D*N) div al) eq 1;
	if (al ne 1) then
	    total_e +:= NumFixedPoints(D, N, al);
	end if;
    end for;
    if #als eq 1 then 
	s := 0; 
    else
	is_prime_power, two, s := IsPrimePower(#als);
	assert is_prime_power and (two eq 2);
    end if;
    g_big := GenusShimuraCurve(D, N);
    g := 1 + (g_big - 1)/2^s - total_e/2^(s+1);
    assert IsIntegral(g);
    return Floor(g);
end function;

function exp_to_Q(e, N, ps)
    ZZ := Integers();
    return &*[ZZ | ps[i]^Valuation(N,ps[i]) : i in [1..#ps] | e[i] eq 1];
end function;

function ALSubgroups(N)
    ZZ := Integers();
    Qs_in_grp := AssociativeArray();
    ps := PrimeDivisors(N);
    W := AbelianGroup([2 : i in [1..#ps]]);
    subs_W := Subgroups(W);
    Qs := {};
    for H in subs_W do
	exps := {Eltseq(W!h) : h in H`subgroup};
	Include(~Qs, {exp_to_Q(e,N,ps) : e in exps});
    end for;
    return Qs;
end function;

function code_we_ran()
    pairs := FindPairs();
    // relevant_pairs := [p : p in pairs | IsEven(Omega(p[1]))];
    relevant_pairs := pairs;
    genera := [GenusShimuraCurve(a[1], a[2]) : a in relevant_pairs];
    dont_know_idxs := [i : i in [1..#genera] | genera[i] ge 3];
    dont_know_pairs := [relevant_pairs[i] : i in dont_know_idxs];
    // Minimal number of Atkin-Lehners we should quotient by to get something 
    // where there is a chance of hyperelliptic
    al_nums := [MinimalNumberOfALinQuotient(a[1], a[2]) : a in relevant_pairs];
    al_nums_shimura := [MinimalNumberOfALinQuotient(a[1], a[2]) : a in relevant_pairs | a[1] gt 1];
    num_curves := &+[#Divisors(a[1]*a[2]) : a in relevant_pairs];
    D := 6;
    // N := 40; // complicated JL
    N := 55;
    min_num := al_nums[Index(pairs, <D,N>)];
    al_subs := ALSubgroups(D*N);
    al_subs_allowed := [S : S in al_subs | #S ge 2^min_num];
    genera := [];
    for als in al_subs_allowed do
	g := GenusShimuraCurveQuotient(D, N, als);
	Append(~genera, <D, N, als, g>); 
    end for;
end function;
		      
function GetGenera(pairs)
    // Restrict first to coprime D and N
    pairs := [p : p in pairs | GCD(p[1], p[2]) eq 1];
    genera := [];
    for i->p in pairs do
	D := p[1];
	N := p[2];
	min_num := al_nums[i];
	al_subs := ALSubgroups(D*N);
	al_subs_allowed := [S : S in al_subs | #S ge 2^min_num];
	for als in al_subs_allowed do
	    g := GenusShimuraCurveQuotient(D, N, als);
	    Append(~genera, <D, N, als, g>); 
	end for;
	print "i = ", i;
    end for;
    return genera;
end function;

/*
function sum_n_powers(a_p, p, n)
    if n eq 0 then 
	return 2;
    end if;
    if n eq 1 then 
	return a_p;
    end if;
    return a_p*sum_n_powers(a_p, p, n-1) - p*sum_n_powers(a_p, p, n-2);
end function;
*/

function sum_n_powers(mfs, p, n, BV)
    assert n ge 1;
    assert Level(mfs) mod p ne 0;
    T_p_n := HeckeOperator(mfs, p^n);
    T_p_n := ChangeRing(T_p_n, Rationals());
    T_p_n := Solution(BV, BV*T_p_n);
    if n eq 1 then 
	return Trace(T_p_n);
    end if;
    T_p_n_2 := HeckeOperator(mfs, p^(n-2));
    T_p_n_2 := ChangeRing(T_p_n_2, Rationals());
    T_p_n_2 := Solution(BV, BV*T_p_n_2);
    return Trace(T_p_n - p*T_p_n_2);
end function;

// Returns false if X is not subhyperelliptic
// If returns true we don't know (compare point counts)
				
function checkSubhyperelliptic(X)
    g := X[4];
    assert g ge 3;
    N := X[2];
    assert X[1] eq 1;
    ws := [w : w in X[3] | w ne 1];
//    assert X[3] eq {1};
    mfs := CuspForms(N);
    als := [AtkinLehnerOperator(mfs,w) : w in ws];
    V := VectorSpace(Rationals(), Dimension(mfs));
    for al in als do
	V := V meet Kernel(al-1);
    end for;
    BV := BasisMatrix(V);
    ps := [p : p in PrimesUpTo(4*g^2) | N mod p ne 0];
    for p in ps do
	v := 1;
	while (p^v le 4*g^2) do
	    num_pts := p^v+1 - sum_n_powers(mfs, p, v, BV);
	    if (num_pts gt 2*(1+p^v)) then
		// print p, v;
		return false;
	    end if;
	    v +:= 1;
	end while;
    end for;
    return true;
end function;



// [HH] Hasegawa, Hashimoto, "Hyperelliptic modular curves X_0^*(N) 
// with square-free levels
//
// [Ogg] Real points on Shimura Curves
