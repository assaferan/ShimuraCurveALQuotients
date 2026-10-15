// From expansions to an equation: the double cover of the star curve.  Load after vw3.m and vwal.m.
//
// THE SHAPE OF THE ANSWER.  Every model this project stores for an index-2 Atkin-Lehner subgroup
// W of the full group W* is a double cover  y^2 = f(t)  of the star curve X* = X_0^D(N)/W*,
// written over a hauptmodul t of X*.  So two functions have to come out of the expansions: the
// hauptmodul itself, and a function on X/W that changes sign under the deck transformation of
// X/W -> X*.  Both are ratios of forms of the SAME weight, so the (1-w)^k automorphy factors
// cancel and only the bare power series is needed.
//
//   t = f1/f2    with f1, f2 a basis of a TWO-dimensional space of forms on the star;
//   y = g/h      with g anti-invariant under the deck element and h invariant, same weight.
//
// y^2 is invariant under the whole of W*, hence a rational function of t, and fitting
// P(t) - y^2 Q(t) = 0 at sample points recovers it.  The curve is then Y^2 = P(t)Q(t) with
// Y = y Q(t), and the square factors of P*Q are dropped.
//
// WHAT THIS CANNOT SEE.  The construction is complex-analytic, so it determines the equation only
// up to a change of coordinate on X* and a scalar on y -- that is, up to isomorphism over C.  The
// QUADRATIC TWIST is invisible to it: y and i*y give the same curve over C and f and -f do not
// give the same curve over Q.  The twist has to come from counting points, which is arithmetic.
//
// THE ANTI-INVARIANT ROWS.  The Gamma and Atkin-Lehner rows of vw3.m are  w'^n - J^k w^n,
// expressing f|g = f.  Asking instead for f|g = -f gives  w'^n + J^k w^n.  So one sign per block
// of samples is all the extra machinery the eigenspace splitting needs.

// Blocks are <samples, eps> with eps = 1 for an invariance relation and -1 for an anti-invariance
// one; samples are the <w, w', matrix> triples that ReducedSamples and ALSamples produce.
function SystemFrom(blocks, k, nterms, CC)
    rows := [];
    for B in blocks do
        eps := B[2];
        for s in B[1] do
            w := s[1]; wp := s[2]; jk := eps * DiscJ(s[3], w)^k;
            Append(~rows, [ wp^n - jk*w^n : n in [0..nterms-1] ]);
        end for;
    end for;
    return Matrix(CC, #rows, nterms, &cat rows);
end function;

// Rows from samples that were reduced over the NORMALISER, each carrying the Hall divisor
// labelling the coset its reduction landed in.  A form invariant under the Atkin-Lehner subgroup
// W, with eigenvalue eps[m] under the rest, satisfies f|g = eps[label] f, so one set of samples
// serves the invariant space, the star space and the anti-invariant space at once: only the sign
// table changes.  A label with no entry in the table contributes no row, which is how the plain
// W-invariant space is cut out from samples that also visited the other coset.
// ⚠ EACH ROW IS NORMALISED TO UNIT LENGTH, and that is not cosmetic.  Reducing over the
// normaliser produces composite group elements whose automorphy factor J(g,w) is far from 1 --
// at X_0^35(1) the weight-4 factor ranges over eight orders of magnitude across the samples --
// and a row scales with J^k.  Equilibrating the COLUMNS, which is what the series expansion
// needs, does nothing about that, so without this the few largest rows dominate every column and
// the system comes out rank deficient: at X_0^35(1) the star space produced a smooth ramp with
// its smallest singular values twenty orders of magnitude BELOW the truncation floor, which is
// the ill-conditioning signature and not a resolution failure.  Scaling a row does not change
// the kernel, and it also makes each relation carry the same RELATIVE truncation error.
// ⚠ WHEN THERE ARE MORE ELIGIBLE SAMPLES THAN ROWS WANTED, THEY ARE TAKEN WITH A STRIDE AND NOT
// AS A PREFIX.  The samples come in order of angle around the sampling circle, so keeping the
// first nmax of them keeps a contiguous ARC -- a sixth of the circle at the row budgets used
// here -- and fitting a polynomial of degree NT on a short arc is exponentially ill-conditioned
// however many points it has.  That is the same defect, in a new place, as the one that made the
// larger bases look unreachable: there the arcs came from a reduction that silently failed, here
// from a cap written as a break.  The symptom is identical and worth recognising: a smooth
// geometric ramp of singular values running far BELOW the truncation floor, with no gap.
//
// Each row is also normalised to unit length, so every relation carries the same relative
// truncation error whatever its automorphy factor.
function SystemLabelled(samples, k, nterms, sgn, nmax, CC)
    elig := [ s : s in samples | IsDefined(sgn, s[4]) ];
    step := Max(1, #elig div Max(1, nmax));
    rows := []; jlo := -1.0; jhi := 0.0;
    i := 1;
    while i le #elig and #rows lt nmax do
        s := elig[i]; i +:= step;
        w := s[1]; wp := s[2]; j := Abs(DiscJ(s[3], w))^k;
        if jlo lt 0 or j lt jlo then jlo := j; end if;
        if j gt jhi then jhi := j; end if;
        jk := sgn[s[4]] * DiscJ(s[3], w)^k;
        r := [ wp^n - jk*w^n : n in [0..nterms-1] ];
        nm := Sqrt(&+[ Abs(x)^2 : x in r ]);
        if nm le 0 then continue; end if;
        Append(~rows, [ x/nm : x in r ]);
    end while;
    return Matrix(CC, #rows, nterms, &cat rows), #rows, jlo, jhi;
end function;

// Rows asking a form to vanish at given points of the disc.  A space of weight-k forms on a
// GENUS-0 quotient is the sections of a line bundle of degree dim - 1 on the line, so a ratio of
// two of them is a map of degree dim - 1, and a hauptmodul needs dim - 2 common zeros imposed.
// That is what these rows do: each point of the disc is one linear condition, the same condition
// at every point of its orbit, so dim - 2 generic points leave a pencil whose members share
// exactly those zeros and whose ratio has degree one.  The (1-w)^k factor never vanishes in the
// disc, so vanishing of the series is vanishing of the form.
function VanishRows(pts, nterms, CC)
    rows := [];
    for w in pts do
        r := [ CC | w^n : n in [0..nterms-1] ];
        nm := Sqrt(&+[ Abs(x)^2 : x in r ]);
        Append(~rows, [ x/nm : x in r ]);
    end for;
    return Matrix(CC, #rows, nterms, &cat rows);
end function;

// The d smallest right singular vectors, with the column equilibration undone, together with the
// spectrum.  A complex singular value appears TWICE in the real form of the matrix, so the vectors
// are taken from every second row from the bottom.  Returns the basis, the ascending singular
// values, and the gap between the d-th and the (d+1)-th -- the number that says whether the space
// was resolved at all.
function KernelBasis(Mx, d, CC)
    n := Ncols(Mx);
    sc := [ Sqrt(&+[ Abs(Mx[i][j])^2 : i in [1..Nrows(Mx)] ]) : j in [1..n] ];
    for j in [1..n] do
        if sc[j] gt 0 then for i in [1..Nrows(Mx)] do Mx[i][j] := Mx[i][j]/sc[j]; end for; end if;
    end for;
    MR := RealForm(Mx);
    Sig, L, Rt := SingularValueDecomposition(MR);
    r := Min(Nrows(MR), Ncols(MR));
    svall := Reverse([ Abs(Sig[i][i]) : i in [1..r] ]);
    sv := [ svall[i] : i in [1..#svall by 2] ];
    bas := [];
    for a in [0..d-1] do
        row := r - 2*a;
        Append(~bas, [ CC | (Rt[row][i] + CC.1*Rt[row][n+i])/sc[i] : i in [1..n] ]);
    end for;
    gap := (d lt #sv and sv[d] gt 0) select sv[d+1]/sv[d] else 0.0;
    return bas, sv, gap;
end function;

// sum b_n w^n -- the (1-w)^k factor is omitted because every use here is a ratio of two forms of
// the same weight, where it cancels.
function Ser(bs, w)
    s := Parent(w) ! 0; p := Parent(w) ! 1;
    for n in [1..#bs] do s +:= bs[n]*p; p *:= w; end for;
    return s;
end function;

// Fit  P(t) = y^2 Q(t)  with deg P <= dp and deg Q <= dq, from values of t and y at sample points.
// Returns the coefficient vectors and the two smallest singular values; the fit is believable when
// the first is at the truncation level and the second is of order 1.
function FitRational(ts, ys, dp, dq, CC)
    rows := [];
    for i in [1..#ts] do
        t := ts[i]; y2 := ys[i]^2;
        r := [ CC | t^j : j in [0..dp] ] cat [ CC | -y2*t^j : j in [0..dq] ];
        nm := Sqrt(&+[ Abs(x)^2 : x in r ]);
        Append(~rows, [ x/nm : x in r ]);
    end for;
    M := Matrix(CC, #rows, dp+dq+2, &cat rows);
    v, s := SmallestVec(M);
    sv := SingVals(M);
    return [ v[j] : j in [1..dp+1] ], [ v[dp+1+j] : j in [1..dq+1] ], sv[1], sv[2];
end function;

// Roots of a complex polynomial given by its coefficient list (ascending), with roots that agree
// to within tol merged and their multiplicity recorded.  Used to drop the square factors of P*Q:
// the branch points are the roots of ODD multiplicity.
function ClusteredRoots(cs, tol, CC)
    while #cs gt 1 and Abs(cs[#cs]) lt tol * Max([Abs(c) : c in cs]) do Prune(~cs); end while;
    PC := PolynomialRing(CC);
    f := &+[ cs[i]*PC.1^(i-1) : i in [1..#cs] ];
    rs := [ z[1] : z in Roots(f) ];
    cls := [];
    for z in rs do
        done := false;
        for i in [1..#cls] do
            if Abs(cls[i][1] - z) lt tol then
                cls[i] := <(cls[i][1]*cls[i][2] + z)/(cls[i][2]+1), cls[i][2]+1>; done := true; break;
            end if;
        end for;
        if not done then Append(~cls, <z, 1>); end if;
    end for;
    return cls, LeadingCoefficient(f);
end function;

// Does ONE Mobius map in t carry a published set of branch points onto the computed one?  The
// comparison a complex-analytic construction can make against a published equation: the two
// coordinates on X* differ by a Mobius change, so the branch sets must agree as sets up to one
// such change.  Three of the published roots fix the map and the rest are predictions, so the
// orderings are tried exhaustively (at most 6*5*4 of them) and the best residual reported.
function MobiusMatch(comp, pub, CC)
    n := #comp;
    if n ne #pub or n lt 4 then return false, 0.0, []; end if;
    best := Infinity(); bestperm := [];
    for i in [1..n] do for j in [1..n] do for l in [1..n] do
        if #{i,j,l} ne 3 then continue; end if;
        // the Mobius map sending pub[i],pub[j],pub[l] to comp[1],comp[2],comp[3].  The normalised
        // coordinate sends the third reference point to infinity, so a point equal to it is
        // reported as a large number rather than dividing by zero.
        cr := func< z, a, b, c | Abs((z-c)*(b-a)) lt 1e-30 select Parent(z) ! 10^30
                                 else ((z-a)*(b-c))/((z-c)*(b-a)) >;
        ok := true;
        img := [];
        for m in [1..n] do
            // match by cross-ratio: send everything to the normalised coordinate
            u := cr(pub[m], pub[i], pub[j], pub[l]);
            Append(~img, u);
        end for;
        tgt := [ cr(comp[m], comp[1], comp[2], comp[3]) : m in [1..n] ];
        // every normalised published value must be matched by a normalised computed one
        tot := 0.0; used := {};
        for u in img do
            dmin := Infinity(); km := 0;
            for m in [1..n] do
                if m in used then continue; end if;
                dd := Abs(u - tgt[m])/(1 + Abs(u));
                if dd lt dmin then dmin := dd; km := m; end if;
            end for;
            if km eq 0 then ok := false; break; end if;
            Include(~used, km); tot +:= dmin;
        end for;
        if ok and tot lt best then best := tot; bestperm := [i,j,l]; end if;
    end for; end for; end for;
    return best lt 1e-6, best, bestperm;
end function;
