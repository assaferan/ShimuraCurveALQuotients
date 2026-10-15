// Molina's algorithm for equations of hyperelliptic Shimura curves of odd discriminant and of
// their Atkin-Lehner quotients.   S. Molina, "Equations of hyperelliptic Shimura curves",
// arXiv:1004.3675v1, Sections 7 and 10.
//
// Load after AttachSpec("ShimuraQuotients.spec").
//
// The model sought is  W : y^2 = R(x) = prod_i p_{R_i}(x)  over Z[1/2], deg R = 2g+2, where the
// Weierstrass points of W are the Heegner points fixed by the involution that becomes hyperelliptic,
// grouped into Galois orbits CM(R_i).  The four steps supply, respectively: the degree of each
// p_{R_i} and the valuations of its discriminant; the leading coefficients; the number field each
// p_{R_i} defines; and finally the polynomials themselves, from IndexFormEquation.

Zx<x> := PolynomialRing(Integers());
Qx<xQ> := PolynomialRing(Rationals());

// ---------------------------------------------------------------- Step 1 (combinatorial part)

// Prop 3.1: the quadratic orders whose Heegner points are fixed by w_k, by discriminant.
function MolFixedDiscs(k)
    if k eq 1 then return []; end if;
    if k eq 2 then return [-4, -8]; end if;
    if (k mod 4) eq 3 then return [-k, -4*k]; end if;
    return [-4*k];
end function;

// Y = X_0^D / <w_m> (m = 1: the curve itself), assumed hyperelliptic of genus g.  Returns the
// divisors n | D whose involution induces the hyperelliptic involution of Y, each with the
// factorisation data of R(x): a list of <disc of R_i, degree of p_{R_i}, number of such factors>.
// Section 10.2: WP(Y) is the image of Fix(w_n) u Fix(w_n'), n' = n m / gcd(n,m)^2, and the degree
// of each factor is [Q(pi(P)) : Q], Prop 10.1 -- which is Gonzalez-Rotger Thm 5.12, the repo's
// DegreeOfFieldOfDefinitionOfCMPoint.
function MolWeierstrassData(D, m, g)
    W := (m eq 1) select {Integers() | 1} else {Integers() | 1, m};
    X := CreateShimuraQuot(D, 1, W);
    out := [];
    for n in Divisors(D) do
        if n eq m or n eq 1 then continue; end if;
        np := (n*m) div GCD(n, m)^2;
        if np lt n then continue; end if;
        discs := Sort(SetToSequence(Seqset(MolFixedDiscs(n) cat MolFixedDiscs(np))));
        tot := 0; rows := [];
        for d in discs do
            nemb := NumberOfOptimalEmbeddings(QuadraticOrder(BinaryQuadraticForms(d)), D, 1);
            if nemb eq 0 then continue; end if;
            deg := DegreeOfFieldOfDefinitionOfCMPoint(X, d);
            if deg eq 0 then continue; end if;
            npts := (m eq 1 or d in MolFixedDiscs(m)) select nemb else nemb div 2;
            if npts mod deg ne 0 then continue; end if;
            tot +:= npts;
            Append(~rows, <d, deg, npts div deg>);
        end for;
        if tot eq 2*g + 2 then Append(~out, <n, np, rows>); end if;
    end for;
    return out;
end function;

// ---------------------------------------------------------------- Step 1 (discriminant of R)

// Thicknesses of the singular points of the special fibre at p | D of X_0^D / W.
// The edges of the dual graph are the singular points of X_p, of thickness #(O^x / {+-1})
// (Lemma 4.1).  Section 10.1: a w_m with p | m swaps the two branches at an edge it fixes, so the
// image is a smooth point; a w_m with p not dividing m fixes both branches and the thickness
// doubles; edges moved by W pair up and keep their thickness.
function MolQuotientThicknesses(D, p, W)
    data := DualGraphData(D, 1, p);
    nE := data[2]; alE := data[4]; uE := data[7];
    act := AssociativeArray();
    act[1] := [1 .. nE];
    for t in alE do act[t[1]] := t[2]; end for;
    for mm in Sort(Setseq(W)) do
        if IsDefined(act, mm) then continue; end if;
        img := [1 .. nE];
        for l in PrimeDivisors(mm) do
            q := l^Valuation(D, l);
            img := [act[q][img[e]] : e in [1 .. nE]];
        end for;
        act[mm] := img;
    end for;
    out := []; seen := {};
    for e in [1 .. nE] do
        if e in seen then continue; end if;
        seen join:= {act[mm][e] : mm in W};
        stab := {mm : mm in W | act[mm][e] eq e};
        if exists{mm : mm in stab | mm mod p eq 0} then continue; end if;
        Append(~out, (uE[e] div 2) * #stab);
    end for;
    return out;
end function;

// |disc R(x)|.  At p | D it is the sum of the thicknesses over the singular points of the fibre
// (Thm 2.3 summed over the pairs of Weierstrass points that collide); at 2 it is 4(g+1), since
// X_0^D has good reduction at 2 for odd D (equation (7.30)); nowhere else, because the model is
// smooth away from 2D.
function MolDiscR(D, m, g)
    W := (m eq 1) select {Integers() | 1} else {Integers() | 1, m};
    v := AssociativeArray();
    th := AssociativeArray();
    v[2] := 4*(g+1);
    th[2] := [];
    for p in PrimeDivisors(D) do
        t := Reverse(Sort(MolQuotientThicknesses(D, p, W)));
        th[p] := t;
        v[p] := &+t;
    end for;
    return v, &*[Integers() | p^v[p] : p in Keys(v)], th;
end function;

// Upper bound for nu_p(disc p_i) at p | D: each singular point of the fibre carries exactly two
// Weierstrass points (Thm 2.3), so a factor of degree n can account for at most floor(n/2) of
// them, and contributes at most the sum of the floor(n/2) largest thicknesses.
function MolDiscBound(th, p, n, vp)
    if p eq 2 then return vp; end if;
    k := n div 2;
    t := th[p];
    if k eq 0 then return 0; end if;
    return &+t[1 .. Min(k, #t)];
end function;

// ---------------------------------------------------------------- Step 2 (leading coefficient)

// Is N the norm of an element of the maximal order of discriminant d (class number 1)?
function MolIsNorm(N, d)
    if N lt 0 then return false; end if;
    if N eq 0 then return true; end if;
    return &and[IsEven(t[2]) : t in Factorization(N) | KroneckerSymbol(d, t[1]) eq -1];
end function;

// The intersection index (P, w_m(P))_p of (6.25), for P in CM(R_infty) with R_infty the maximal
// order of discriminant d and class number 1.  The norm equation (6.26) is
//     (-d) m = n + N(lambda_-) D p,     n = (-d) N(lambda_+) >= 0,
// with n and N(lambda_-) norms from K_infty; the index is then
//     max ord_p(N(lambda_-)) / (1 - (d|p)) + 1,
// and 0 when there is no solution.  Only primes p that do not split in K_infty can occur, since
// P_infty has supersingular specialisation at every prime dividing the leading coefficient.
function MolIntersectionIndex(D, d, m, p)
    if D mod p eq 0 then return 0; end if;
    chi := KroneckerSymbol(d, p);
    if chi eq 1 then return 0; end if;
    best := -1;
    for Nm in [1 .. ((-d)*m) div (D*p)] do
        n := (-d)*m - Nm*D*p;
        if not (MolIsNorm(n, d) and MolIsNorm(Nm, d)) then continue; end if;
        val := Valuation(Nm, p);
        if val mod (1 - chi) eq 0 then best := Max(best, val div (1 - chi)); end if;
    end for;
    return best lt 0 select 0 else best + 1;
end function;

// The imaginary quadratic fields usable for the points at infinity: maximal, class number 1,
// discriminant prime to D, and with CM(R_infty) non-empty (Section 3, Section 10.2).
function MolInfinityDiscs(D)
    return [d : d in [-3,-4,-7,-8,-11,-19,-43,-67,-163] |
            GCD(d, D) eq 1 and
            NumberOfOptimalEmbeddings(QuadraticOrder(BinaryQuadraticForms(d)), D, 1) gt 0];
end function;

// The leading coefficient a_R of R(x), for the points at infinity {pi(P_inf), w(pi(P_inf))} with
// P_inf a Heegner point by the maximal order of discriminant d.
// Sign and square class: Q(sqrt(a_R)) is the field of definition of pi(P_inf).
// Valuations: nu_p(a_R) = (1 - (d|p)) * ( (P,w_n P)_p + (P,w_n' P)_p ) by (10.32), the second term
// being absent on the curve itself.
function MolLeadingCoefficient(D, m, n, np, d)
    W := (m eq 1) select {Integers() | 1} else {Integers() | 1, m};
    Y := CreateShimuraQuot(D, 1, W);
    F := FieldOfDefinitionOfCMPoint(Y, d);
    if Degree(F) eq 1 then
        sq := 1;
    else
        sq := Squarefree(Integers() ! Discriminant(MaximalOrder(AbsoluteField(F))));
    end if;
    ms := (m eq 1) select [n] else [n, np];
    a := 1;
    for p in PrimesUpTo(Abs(d)) do
        if D mod p eq 0 then continue; end if;
        chi := KroneckerSymbol(d, p);
        if chi eq 1 then continue; end if;
        idx := &+[MolIntersectionIndex(D, d, mm, p) : mm in ms];
        if idx eq 0 then continue; end if;
        a *:= p^((1 - chi)*idx);
    end for;
    if sq lt 0 then a := -a; end if;
    ok := IsSquare(a / sq);
    return a, ok, sq;
end function;

// ---------------------------------------------------------------- Step 3/4 (the polynomials)

// The IndexFormEquation index attached to a factor of degree n, leading coefficient a and
// discriminant J^2 disc(M).  With q(x) = a^(n-1) p(x/a) monic, disc(q) = a^((n-1)(n-2)) disc(p),
// so the index of its root in O_M is |a|^((n-1)(n-2)/2) J.
function MolIndex(a, n, J)
    return Abs(a)^(((n-1)*(n-2)) div 2) * J;
end function;

// All degree-n integer polynomials with leading coefficient a, defining the field M, and with
// discriminant J^2 disc(M), up to translation by an integer.  Step 4: q(x) = a^(n-1) p(x/a) is
// monic with a root in O_M of the index above, which IndexFormEquation enumerates.
function MolFactorCandidates(M, O, a, n, J : Verbose := false)
    if n eq 1 then
        return [Zx | a*x + b : b in [0 .. Abs(a)-1]];
    end if;
    I := MolIndex(a, n, J);
    if n eq 2 then
        // a generator of index I is u + I*omega for omega a generator of O_M
        w := O.2;
        gens := [M ! (I*w)];
    else
        gens := [M ! s : s in IndexFormEquation(O, I)];
    end if;
    // IndexFormEquation returns its generators only up to sign and translation by Z, but
    // p(x) -> p(x+r) moves the root of q by a*r, so the translates must be put back by hand:
    // only one residue class modulo a can give a p with integer coefficients.
    out := [];
    for gm in gens do
        for sg in [1, -1] do
            for t in [0 .. Abs(a) - 1] do
                g := sg*gm + t;
                q := MinimalPolynomial(g);
                if Degree(q) ne n then continue; end if;
                c := Coefficients(q);                   // c[j+1] is the coefficient of x^j
                good := true;
                cf := [];
                for j in [0 .. n] do
                    e := n - 1 - j;
                    cj := Integers() ! c[j+1];
                    if e le 0 then
                        Append(~cf, cj * a^(-e));
                    else
                        if cj mod Abs(a)^e ne 0 then good := false; break; end if;
                        Append(~cf, ExactQuotient(cj, a^e));
                    end if;
                end for;
                if not good then continue; end if;
                p := Zx ! cf;
                assert LeadingCoefficient(p) eq a;
                Append(~out, p);
            end for;
        end for;
    end for;
    return Sort(SetToSequence(Seqset(out)));
end function;

// Integers r with |Res(f(x), h(x+r))| at most bnd.  All complex roots of f and h lie in a disc of
// radius C about 0, so |Res| >= |lc(f)^deg h lc(h)^deg f| (|r| - 2C)^(deg f deg h).
function MolTranslationRange(f, h, bnd)
    rts := [Abs(z[1]) : z in Roots(f, ComplexField(30))] cat
           [Abs(z[1]) : z in Roots(h, ComplexField(30))];
    C := #rts eq 0 select 0 else Ceiling(Maximum(rts));
    lc := Abs(LeadingCoefficient(f))^Degree(h) * Abs(LeadingCoefficient(h))^Degree(f);
    e := Degree(f) * Degree(h);
    R := 2*C + Ceiling((bnd / lc)^(1/e)) + 2;
    return [-R .. R];
end function;
