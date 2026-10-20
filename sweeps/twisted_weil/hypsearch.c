/* Exhaustive search for genus-3 hyperelliptic curves y^2 = f(x) over F_q (q odd prime) with given (N1, N2).
   Every such curve has a model with deg f = 8 (P^1(F_q) has q+1 > 8 points, so move a non-branch point to oo),
   and then (shift x, scale x, scale y) f = c (x^8 + a6 x^6 + a5 x^5 + ... + a0) with c in {1, n}, a6 in {0, 1, n},
   n a fixed non-residue.  Prints every (c, a6..a0) with #C(F_q) = N1 and #C(F_q^2) = N2 (no squarefree check:
   singular f are filtered afterwards in Magma, which also computes the exact Weil polynomial).
   usage: hypsearch q N1 N2 */
#include <stdio.h>
#include <stdlib.h>
int main(int argc, char **argv) {
  int q = atoi(argv[1]); long N1t = atol(argv[2]), N2t = atol(argv[3]);
  int chi[64]; for (int i = 0; i < q; i++) chi[i] = -1; chi[0] = 0;
  for (int i = 1; i < q; i++) chi[(i * i) % q] = 1;
  int n = 0; for (int i = 2; i < q; i++) if (chi[i] == -1) { n = i; break; }
  int cs[2] = {1, n}, a6s[3] = {0, 1, n};
  long found = 0, n1pass = 0;
  int a[9];
  for (int ci = 0; ci < 2; ci++) for (int ai = 0; ai < 3; ai++) {
    int c = cs[ci];
    a[8] = 1; a[7] = 0; a[6] = a6s[ai];
    for (long idx = 0; idx < (long)q*q*q*q*q*q; idx++) {
      long t = idx; for (int k = 0; k <= 5; k++) { a[k] = t % q; t /= q; }
      /* N1 */
      long N1 = 1 + chi[c];
      for (int x = 0; x < q; x++) {
        int v = 1; for (int k = 7; k >= 0; k--) v = (v * x + a[k]) % q;
        v = (v * c) % q; N1 += 1 + chi[v];
      }
      if (N1 != N1t) continue;
      n1pass++;
      /* N2 over F_q^2 = F_q[s]/(s^2 - n); z = u + w s; chi2(z) = chi(norm z), norm = u^2 - n w^2 */
      long N2 = 2;
      for (int u = 0; u < q; u++) for (int w = 0; w < q; w++) {
        int vu = 1, vw = 0;   /* Horner */
        for (int k = 7; k >= 0; k--) {
          int nu = (vu * u + n * ((vw * w) % q) + a[k]) % q;
          int nw = (vu * w + vw * u) % q;
          vu = nu; vw = nw;
        }
        vu = (vu * c) % q; vw = (vw * c) % q;
        int nm = ((vu * vu - n * ((vw * vw) % q)) % q + q) % q;
        N2 += 1 + chi[nm];
      }
      if (N2 != N2t) continue;
      found++;
      printf("C %d", c); for (int k = 0; k <= 8; k++) printf(" %d", a[k]); printf("\n");
    }
  }
  fprintf(stderr, "q=%d n=%d N1=%ld N2=%ld n1pass=%ld found=%ld\n", q, n, N1t, N2t, n1pass, found);
  return 0;
}
