// Loerr and upper bounds for number of points in reduction mod p, 
// scaled by 12/p-1, see [HH]

function LowerBound(D, N)
    ps := PrimeDivisors(N);
    c_D := EulerPhi(D);
    c_N := N * &*[Rationals() | 1 + 1/p : p in ps];
    omega := #ps;
    return c_D * c_N / (2^omega);
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

procedure VerifyBound(: r := 6)
    ps := PrimesUpTo(100);
    assert LowerBound(1, &*ps[1..r]) gt UpperBound(ps[r+1]);
    assert LowerBound(&*ps[1..r], 1) gt UpperBound(ps[r+1]);
end procedure;

function FindMaximalN(r)
    ps := PrimesUpTo(100);
    return Floor(UpperBound(ps[r])*2^(r-1));
end function;

function FindMaximalD(r)
    return 1290;
    ps := PrimesUpTo(100);
    ub := Ceiling(UpperBound(ps[r]));
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

function Omega(n)
    return #PrimeDivisors(n);
end function;

function FindPairs()
    pairs := [];
    N0 := FindMaximalN(6);
    D0 := FindMaximalD(6);
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
    
// Quotient by w_m, m divides DN, following [Ogg]
function GenusShimuraCurveQuotient(D, N, m)
    if (m eq 2) then
	orders := [MaximalOrder(QuadraticField(-1)), 
		   MaximalOrder(QuadraticField(-2))];
    elif (m mod 4 eq 3) then
	O<alpha> := MaximalOrder(QuadraticField(-m));
	orders := [O, sub<O | 1, 2*alpha> ];
    else
	orders := [MaximalOrder(QuadraticField(-m))];
    end if;
    g := 0;
    for R in orders do
	h := PicardNumber(R);
	// use Pete's formula (or Ogg)
    end for;
end function;

function GenusShimuraCurveQuotient(D, N, als)
    
end function;

// [HH] Hasegawa, Hashimoto, "Hyperelliptic modular curves X_0^*(N) 
// with square-free levels
//
// [Ogg] Real points on Shimura Curves
