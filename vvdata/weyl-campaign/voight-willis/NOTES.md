# Voight-Willis at discriminant 6 — WORKS.  2026-10-15

Supersedes the earlier no-go in ../vw/NOTES.md **on the method**.  Every number in that file was
invalid (see "The bug").  Scaling to large levels is a separate question, settled elsewhere.

## The paper

Obtained: J. Voight and J. Willis, "Computing power series expansions of modular forms",
https://jvoight.github.io/articles/ModFormPaper-090317.pdf (author's copy; published in
Computations with Modular Forms, 2014).  Implementation follows Sec. 3.

Method as published: f(z) = (1-w)^k sum_{n>=0} b_n w^n, eq (4), with w = (z-p)/(z-pbar) the Cayley
map centred at the Dirichlet-domain centre p.  rho = max{|w| : w in D}.  Sample w on |w| = rho,
REDUCE each into the domain to get g with w' = gw in D, and impose eq (9):
   (1-w')^k sum b_n w'^n = j(g,z)^k (1-w)^k sum b_n w^n,  j the HALF-PLANE factor cz+d,
kernel eq (11).  N ~ log(eps)/log(rho), eq (8).

Dividing eq (9) by (1-w')^k gives sum b_n (w'^n - jD^k w^n) = 0 with jD the disc automorphy factor.
**Verified numerically that jD^k equals j(g,z)^k (1-w)^k/(1-w')^k to 1e-60** (dbg_b.m), so the disc
form and the paper's form are the same relation; the (1-w)^k factor is a change of basis, not a
correctness issue.

## The bug in the earlier attempt: Magma's SingularValueDecomposition convention

`SingularValueDecomposition(M)` for M of size m x n returns a triple whose shapes are (m x n,
m x m, n x n).  The **singular values are the diagonal of the FIRST return, in DESCENDING order**;
the second is the left orthogonal matrix and the third holds the right singular vectors **as its
ROWS**.  The earlier `vwlib.m` did `_, S, _ := SingularValueDecomposition(MR)` and read the
diagonal of `S` — the orthogonal matrix.  Those entries are bounded by 1 and otherwise arbitrary,
including negative.

Demonstrated on a random 7 x 3 matrix (shape2.m, shape3.m, shape4.m): true singular values from
sqrt of the eigenvalues of M^T M were 34.0418, 13.4420, 10.0274; the diagonal of the first return
reproduced them exactly, while the diagonal of the second read 0.513, 0.156, -0.084, 0.908, ...
Also |U^T U - I| = 1158 for the second return read as the first, and the right singular vectors
satisfy ||M * (row j of third)|| = sv[j] exactly, while the columns do not.

So **every singular value in ../vw/NOTES.md is meaningless** — the "four exact null vectors at
machine zero", the "truncation-level agreement at k = 2 mod 4", the weight-2 failures, the claim
that the domain radius caps resolution.  The apparent k mod 4 pattern was coincidence in where
small entries of an orthogonal matrix happened to fall.

Two further points that are real:
* **Column equilibration is required.**  Columns scale like rho^n, so without it every singular
  value sits near rho^N and no gap is visible (dims6.log: all ~1e-13 at N = 50).  With it the gaps
  are 10^6 to 10^9 (dims6b.log).
* Sampling follows the paper (radius rho, reduced).  I did **not** re-test whether the earlier
  isometric-circle sampling also works once the SVD is read correctly, so the sampling change is
  not shown to have been necessary.

## Validation at D = 6 (X_0^6(1), signature (0; 2,2,3,3), genus 0, rho = 0.543496)

1. **Dimensions** dim M_k = (k-1)(g-1) + sum floor((k/2)(1-1/e_i)), k = 2,4,6,8,10,12,14:
   expected 0,1,1,1,1,3,1 — **all 7 correct** at prec 100, N = 50, 229 samples, gaps 10^6-10^9.
2. **Valence** (not forced by symmetry): deg div(f_k) = k/6 here, so at k = 4 the value 2/3 can only
   be carried by the two order-3 points at 1/3 each.  The computed weight-4 form vanishes at exactly
   those two and at neither order-2 point: relative |f4| = 8.4e-21, 8.7e-21 against 0.503 and 1.000
   at N = 80 (1.4e-13 at N = 50 — converging).
3. **Elkies' configuration.**  From tests/CMFieldsElkies.m: the function field of X_0^6(1) is
   Q(t)(sqrt(-t), sqrt(3(t-1))) with t = 0,1,oo at the star elliptic points of order 2,4,6
   (disc -24,-4,-3).  Setting y^2 = -t, u^2 = 3(t-1) gives the conic u^2 = -3(y^2+1), parametrised
   by m with y = (m^2-1)/(2m), so **t = -((m^2-1)/(2m))^2** and the four elliptic points of
   X_0^6(1) are m = 0, oo (order 3, disc -3) and m = +-i (order 2, disc -4) — cross-ratio -1
   (predict.m).  Built a Hauptmodul as a ratio of two weight-12 forms sharing an imposed zero
   (dim M_12 = 3, so one extra linear condition leaves a 2-dimensional space and the ratio has
   degree 1).  Computed cross-ratio:
       N = 50: -0.999999993526 - 1.902e-8 i   (error 2.0e-8)
       N = 80: -1.00000000000 - 7.73e-15 i    (error 1.1e-14)
   and after normalising three of the points the four land at exactly 0, oo, i, -i.
4. **Elkies' CM value, the sharp check.**  Elkies gives t = 0 at discriminant -24.  Taking mu in O
   with trace 0 and reduced norm 6 (so Z[mu] = Z[sqrt -6], disc -24) and its fixed point, reduced
   into the domain:
       mu = -i + 3j : m = -1.000000000 - 1.85e-17 i,  t = 2.5e-34
       mu = 3j + k  : m = +1.000000000 - 1.21e-17 i,  t = -5.6e-34
   Both predicted values m = -1 and m = +1 appear, t = 0 to 1e-34.
   (Magma's own `CMPoints` is unusable here: it fails inside
   package/Geometry/GrpPSL2/GrpPSL2Shim/reciprocity.m line 29 with a RingClassGroup coercion error.
   Taking the fixed point of mu directly avoids it.)

## Not done

* The Atkin-Lehner layer: all of this is on X_0^D(N), not on a quotient X/W.  Still extra work
  beyond the published method.
* Whether isometric-circle sampling also works with the SVD read correctly.
* Anything about scaling to the levels of the open curves.

## Files
vw3.m (library, correct SVD), dims.m (dimension ladder), haupt.m (Hauptmodul + cross-ratio +
valence), cmchk.m (disc -24 CM check), predict.m (Elkies prediction), dbg_a/b/c/d/e.m,
shape2/3/4.m (SVD convention).
