// Does the preprint's prop:closedcoef scalar a_E(m) reproduce the MEASURED A_m
// that SchoferFormula.m:593-596 records for the Kappa0 log-N defect?
//
//   A_m = -b(m)/4,  b = weight-3/2 Eisenstein coefficient  (SchoferFormula.m:600)
//   -a_E(m) = 12 H(4m) * prod_{p|D} (1-chi_p)/((p-1)(sig_kp - chi_p sig_{kp-1}))
//                      * prod_{p|N} (p-chi_p) p^kp/((p^2-1)(sig_kp - chi_p sig_{kp-1}))
//
// NOTE ON LOAD ORDER: `import` compiles TraceFormula.m as its own package immediately,
// while AttachSpec is lazy.  Unless a ClassNumberData intrinsic is touched FIRST, the
// imported copy cannot resolve ClassNumberLU and H() dies at TraceFormula.m:192.
AttachSpec("ShimuraQuotients.spec");
_ := ClassNumberLU(-4);                       // force ClassNumberData.m to load
import "TraceFormula.m" : Hurwitz, H;         // the repo's own Hurwitz class number

function sig(j, p)
    if j lt 0 then return 0; end if;
    return &+[p^i : i in [0..j]];
end function;

// hurw = 1 -> use the repo's Hurwitz(); 2 -> use the repo's Kronecker-Hurwitz H()
function aE(m, D, N, hurw)
    D0 := FundamentalDiscriminant(-4*m);
    ok, c := IsSquare((-4*m) div D0);
    if not ok then return 0, "BADC"; end if;
    val := 12 * ((hurw eq 1) select Hurwitz(4*m) else H(4*m));
    for p in PrimeDivisors(D) do
        kp := Valuation(c, p); chi := KroneckerSymbol(D0, p);
        den := (p-1)*(sig(kp,p) - chi*sig(kp-1,p));
        if den eq 0 then return 0, "DIV0"; end if;
        val *:= (1 - chi)/den;
    end for;
    if N gt 1 then
        for p in PrimeDivisors(N) do
            kp := Valuation(c, p); chi := KroneckerSymbol(D0, p);
            den := (p^2-1)*(sig(kp,p) - chi*sig(kp-1,p));
            if den eq 0 then return 0, "DIV0"; end if;
            val *:= (p - chi)*p^kp/den;
        end for;
    end if;
    return val, "ok";        // val is -a_E(m)
end function;

cases := [
  <15, 2, [<2,-1>,   <10,1>,   <30,0>]>,
  <6,  5, [<10,3/2>, <15,1/2>, <30,3/2>]>,
  <10, 3, [<3,1/2>,  <12,1/2>, <30,3/2>]>,
  <21, 2, [<2,-2>,   <6,-4>,   <18,2>,  <42,-8>]>
];

for hurw in [1,2] do
  printf "\n===== Hurwitz source: %o =====\n",
      (hurw eq 1) select "repo Hurwitz() [TraceFormula.m:152]" else "repo H() [TraceFormula.m:174]";
  printf "%-8o %-4o %-10o %-10o %-10o %-10o %-10o %o\n",
      "base","m","H(4m)","A_m(meas)","-a_E(m)","-a_E/4","+a_E/4","verdict";
  nm := 0; nt := 0; nz := 0; nzmeas := 0;
  for cs in cases do
    D := cs[1]; N := cs[2];
    for pr in cs[3] do
      m := pr[1]; meas := pr[2]; nt +:= 1;
      hv := (hurw eq 1) select Hurwitz(4*m) else H(4*m);
      v, st := aE(m, D, N, hurw);
      if st ne "ok" then
        printf "%-8o %-4o %-10o %-10o %-10o %-10o %-10o %o\n",
            Sprintf("%o_%o",D,N), m, hv, meas, st, "-", "-", "ERR";
        continue;
      end if;
      c1 := -v/4; c2 := v/4;
      tag := (c1 eq meas) select "MATCH -a_E/4"
             else ((c2 eq meas) select "MATCH +a_E/4" else "no");
      if tag ne "no" then nm +:= 1; end if;
      if v eq 0 then
          nz +:= 1;
          if meas ne 0 then nzmeas +:= 1; end if;
      end if;
      printf "%-8o %-4o %-10o %-10o %-10o %-10o %-10o %o\n",
          Sprintf("%o_%o",D,N), m, hv, meas, v, c1, c2, tag;
    end for;
  end for;
  printf "MATCHED %o of %o\n", nm, nt;
  printf "a_E vanished %o of %o times; of those, %o had a NONZERO measured A_m\n", nz, nt, nzmeas;
end for;

// Cross-check: do the two repo Hurwitz functions ever disagree on the arguments used here?
printf "\n===== Hurwitz() vs H() on all 4m used above =====\n";
ms := [2,10,30,15,3,12,6,18,42];
for m in ms do
    a := Hurwitz(4*m); b := H(4*m);
    printf "  4m=%-5o Hurwitz=%-8o H=%-8o %o\n", 4*m, a, b,
        (a eq b) select "agree" else "DISAGREE";
end for;
quit;
