/* Exhaustive search for genus-g hyperelliptic curves y^2 = f(x) over F_q (q odd prime) with given (N1, N2).
   Every such curve has a model of degree 2g+2 (if some point of P^1(F_q) is not a branch point, move it to oo)
   or of degree 2g+1 (all rational points are branch points: move one to oo).
   even: f = c (x^{2g+2} + a_{2g} x^{2g} + ... + a0), c in {1,n}, a_{2g} in {0,1,n}  (shift kills x^{2g+1}; needs p | 2g+2 false)
   odd:  f = x^{2g+1} + a_{2g-1} x^{2g-1} + ... + a0                               (shift kills x^{2g}; needs p | 2g+1 false;
         leading coefficient normalised to 1 by x -> c x, y -> scaled)
   Prints candidates "C c a0 .. a_deg"; Magma then drops singular f and computes exact Weil polynomials.
   usage: hypsearch2 q g N1 N2 */
#include <stdio.h>
#include <stdlib.h>
int q, n, chi[64];
long N1t, N2t, found = 0, n1pass = 0;
static void test(int *a, int deg, int c) {
  long N1 = (deg % 2 == 0) ? 1 + chi[c] : 1;
  for (int x = 0; x < q; x++) {
    int v = 0; for (int k = deg; k >= 0; k--) v = (v * x + a[k]) % q;
    v = (v * c) % q; N1 += 1 + chi[v];
  }
  if (N1 != N1t) return;
  n1pass++;
  long N2 = (deg % 2 == 0) ? 2 : 1;
  for (int u = 0; u < q; u++) for (int w = 0; w < q; w++) {
    int vu = 0, vw = 0;
    for (int k = deg; k >= 0; k--) {
      int nu = (vu * u + n * ((vw * w) % q) + a[k]) % q;
      int nw = (vu * w + vw * u) % q;
      vu = nu; vw = nw;
    }
    vu = (vu * c) % q; vw = (vw * c) % q;
    int nm = ((vu * vu - n * ((vw * vw) % q)) % q + q) % q;
    N2 += 1 + chi[nm];
  }
  if (N2 != N2t) return;
  found++;
  printf("C %d", c); for (int k = 0; k <= deg; k++) printf(" %d", a[k]); printf("\n");
}
int main(int argc, char **argv) {
  q = atoi(argv[1]); int g = atoi(argv[2]); N1t = atol(argv[3]); N2t = atol(argv[4]);
  for (int i = 0; i < q; i++) chi[i] = -1; chi[0] = 0;
  for (int i = 1; i < q; i++) chi[(i * i) % q] = 1;
  n = 0; for (int i = 2; i < q; i++) if (chi[i] == -1) { n = i; break; }
  int a[32];
  /* even models */
  if ((2*g+2) % q != 0) {
    int deg = 2*g+2, cs[2] = {1, n}, tops[3] = {0, 1, n};
    long tot = 1; for (int k = 0; k < 2*g; k++) tot *= q;
    for (int ci = 0; ci < 2; ci++) for (int ti = 0; ti < 3; ti++)
      for (long idx = 0; idx < tot; idx++) {
        long t = idx; for (int k = 0; k < 2*g; k++) { a[k] = t % q; t /= q; }
        a[2*g] = tops[ti]; a[2*g+1] = 0; a[2*g+2] = 1;
        test(a, deg, cs[ci]);
      }
  } else {  /* p | 2g+2: cannot shift away x^{2g+1}; normalise a_{2g+1} in {0,1} by x -> lambda x instead */
    int deg = 2*g+2, cs[2] = {1, n};
    long tot = 1; for (int k = 0; k < 2*g+1; k++) tot *= q;
    for (int ci = 0; ci < 2; ci++) for (int ti = 0; ti < 2; ti++)
      for (long idx = 0; idx < tot; idx++) {
        long t = idx; for (int k = 0; k < 2*g+1; k++) { a[k] = t % q; t /= q; }
        a[2*g+1] = ti; a[2*g+2] = 1;
        test(a, deg, cs[ci]);
      }
  }
  /* odd models */
  if ((2*g+1) % q != 0) {
    int deg = 2*g+1;
    long tot = 1; for (int k = 0; k < 2*g; k++) tot *= q;
    for (long idx = 0; idx < tot; idx++) {
      long t = idx; for (int k = 0; k < 2*g; k++) { a[k] = t % q; t /= q; }
      a[2*g] = 0; a[2*g+1] = 1;
      test(a, deg, 1);
    }
  } else { fprintf(stderr, "p | 2g+1: odd normalisation invalid\n"); return 1; }
  fprintf(stderr, "q=%d g=%d n=%d N1=%ld N2=%ld n1pass=%ld found=%ld\n", q, g, n, N1t, N2t, n1pass, found);
  return 0;
}
