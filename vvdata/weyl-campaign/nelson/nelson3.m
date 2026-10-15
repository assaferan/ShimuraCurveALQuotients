// Nelson, Theorem 3.1, third version.  Two changes over nelson2.m, both needed at large levels:
//
//  * the lattice points are produced by a ShortVectorsProcess rather than a list, because at
//    D N = 4830 one truncation holds several million of them and the list alone exhausts memory;
//  * SEVERAL eigenforms are evaluated in one enumeration.  Everything except a_f(n) is shared
//    between the forms of a base, so a basis of g differentials at one pair of points costs one
//    enumeration, not g of them.  That is the shape a canonical model needs.
//
// The W_k table depends only on the weight and the cutoff, not on the base, so it is built once
// by WkTable and handed to every NelsonSetupT.

load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/nelson/wk.m";

function WkDirect(k, x, RR)
    eps := RR ! 10^(-(Precision(RR)+5));
    s := RR ! 0; n := 1;
    while n le 400 do
        u := n*x;
        t := (RR ! n)^k * (u*KBessel(k-1, u) - KBessel(k, u));
        s +:= t;
        if n ge 2 and Abs(t) lt eps*(1 + Abs(s)) then break; end if;
        n +:= 1;
    end while;
    return s / 2^(k-1);
end function;

function WkTable(k, cutmax, RR)
    return WkPanels(func< x | WkDirect(k, x, RR) >, RR ! 0.4, RR ! (cutmax+1), RR
                    : ratio := 1.5, deg := Precision(RR) + 14);
end function;

function TraceGram3(bas, A)
    n := #bas;
    return Matrix(Rationals(), n, n, [ Trace(A ! (bas[i]*bas[j])) : i, j in [1..n] ]);
end function;
function NormGram3(bas, A)
    n := #bas;
    return Matrix(Rationals(), n, n,
        [ Trace(A ! (bas[i]*Conjugate(A ! bas[j])))/2 : i, j in [1..n] ]);
end function;
function RdBasis3(T, DN, d)
    n := Nrows(T);
    gens := VerticalJoin(IdentityMatrix(Rationals(), n), (DN div d)*T^(-1));
    den := LCM([ Denominator(x) : x in Eltseq(gens) ]);
    return ChangeRing(BasisMatrix(Lattice(ChangeRing(den*gens, Integers()))), Rationals())/den;
end function;
function SigZ(z, CC)
    x := Real(z); y := Imaginary(z); r := Sqrt(y);
    return Matrix(CC, 2, 2, [ r, x/r, 0, 1/r ]);
end function;
function SigZi(z, CC)
    x := Real(z); y := Imaginary(z); r := Sqrt(y);
    return Matrix(CC, 2, 2, [ 1/r, -x/r, 0, r ]);
end function;

NelRec3 := recformat< D, N, k, A, O, mrep, DN, ds, latb, latQ, prec, Wtab, cutmax >;

function NelsonSetupT(D, N, k, prec, cutmax, Wtab : OO := 0, MREP := 0)
    if Type(OO) eq RngIntElt then
        A := QuaternionAlgebra(D); O := QuaternionOrder(A, N);
    else
        O := OO; A := Algebra(O);
    end if;
    assert Discriminant(O) eq D*N;
    if Type(MREP) eq RngIntElt then
        mrep := FuchsianMatrixRepresentation(A : Precision := prec);
    else
        mrep := MREP;
    end if;
    bas := [ A ! b : b in ZBasis(O) ];
    T := TraceGram3(bas, A);
    DN := D*N;
    latb := AssociativeArray(); latQ := AssociativeArray();
    for d in Divisors(DN) do
        B := RdBasis3(T, DN, d);
        bd := [ &+[ B[i][j]*bas[j] : j in [1..4] ] : i in [1..4] ];
        latb[d] := bd; latQ[d] := NormGram3(bd, A);
    end for;
    return rec< NelRec3 | D := D, N := N, k := k, A := A, O := O, mrep := mrep, DN := DN,
                          ds := Divisors(DN), latb := latb, latQ := latQ, prec := prec,
                          Wtab := Wtab, cutmax := cutmax >;
end function;

// afs is a sequence of sequences:  afs[i][n] = a_f(n) for the i-th eigenform, n = 1..nmax.
// Returns the sequence of values of Theorem 3.1, one per eigenform, plus the lattice-point count
// and the largest n used.
function NelsonEvalMulti(dat, z1, z2, afs, cut, CC)
    k := dat`k; RR := RealField(Precision(CC)); pi := Pi(RR); i_ := CC.1;
    assert IsEven(k) and cut le dat`cutmax;
    nf := #afs; nmax := #afs[1];
    s1i := SigZi(z1, CC); s2 := SigZ(z2, CC);
    C := cut/(4*pi); Pb := 2*C^2;
    tot := [ CC | 0 : i in [1..nf] ];
    nterms := 0; nmaxused := 0;
    SCL := 10^16;
    for d in dat`ds do
        bs := dat`latb[d]; Qn := dat`latQ[d];
        Ms := [ s1i * ChangeRing(dat`mrep(b), CC) * s2 : b in bs ];
        Xs := [ ((M[1][1]+M[2][2]) + i_*(M[1][2]-M[2][1]))/2 : M in Ms ];
        G := Matrix(RR, 4, 4, [ RR ! 0 : j in [1..16] ]);
        for a in [1..4] do for b in [a..4] do
            v := (&+[ Real(Ms[a][r][c])*Real(Ms[b][r][c]) : r, c in [1..2] ])/2;
            G[a][b] := v; G[b][a] := v;
        end for; end for;
        GI := Matrix(Integers(), 4, 4, [ Round(SCL*G[a][b]) : a, b in [1..4] ]);
        for a in [1..4] do for b in [a+1..4] do GI[b][a] := GI[a][b]; end for; end for;
        for a in [1..4] do GI[a][a] +:= 8; end for;
        L := LatticeWithGram(GI);
        cf := (CC ! MoebiusMu(GCD(d, dat`N))) / (d * afs[1][d gt nmax select 1 else d]);
        cfs := [ (CC ! MoebiusMu(GCD(d, dat`N))) / (d * afs[i][d]) : i in [1..nf] ];
        sub := [ CC | 0 : i in [1..nf] ];
        P := ShortVectorsProcess(L, Ceiling(SCL*Pb*1.000001));
        while not IsEmpty(P) do
            vv := NextVector(P);
            v := Eltseq(vv);
            if v eq [0,0,0,0] then continue; end if;
            t := &+[ v[a]*v[b]*Qn[a][b] : a, b in [1..4] ];
            if t le 0 then continue; end if;
            n := d*t;
            n := Integers() ! n;
            if n gt nmax then continue; end if;
            X := v[1]*Xs[1] + v[2]*Xs[2] + v[3]*Xs[3] + v[4]*Xs[4];
            aX := Abs(X); x := 4*pi*aX;
            if x gt cut then continue; end if;
            wv := (X/aX)^k * WkFast(dat`Wtab, x);
            for i in [1..nf] do sub[i] +:= afs[i][n]*wv; end for;
            nterms +:= 1;
            if n gt nmaxused then nmaxused := n; end if;
        end while;
        for i in [1..nf] do tot[i] +:= 2*cfs[i]*sub[i]; end for;
    end for;
    return tot, nterms, nmaxused;
end function;
