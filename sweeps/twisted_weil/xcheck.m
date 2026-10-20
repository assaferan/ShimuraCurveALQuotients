// Independent recomputation (Eichler-Selberg, TraceDNewALFixed; no modular symbols) of P_{X_h} for AL h = w_Q.
// IN lines: D N g W(comma) Q p expected-key.  Tr(T_n w_Q | S^W) = 2 Tr(T_n | S^<W,Q>) - Tr(T_n | S^W).
// Tr((hF)^k | H^1) = Tr(s_k(T_p) h^k), s_k = T_{p^k} - p T_{p^{k-2}}; #X_h(F_{p^k}) = p^k + 1 - that; Newton -> P.
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
for line in Split(Read(IN), "\n") do
  f := Split(line, " "); if #f lt 7 then continue; end if;
  D, N, g := Explode([StringToInteger(f[i]) : i in [1..3]]);
  W := {StringToInteger(x) : x in Split(f[4], ",")}; Q := StringToInteger(f[5]); p := StringToInteger(f[6]);
  L := D*N; W2 := W join {AtkinLehnerMul(Q, w, L) : w in W};
  tr := function(n, WW) if n lt 1 then return 0; end if; return TraceDNewALFixed(D, N, 2, n, WW); end function;
  S := [];
  for k in [1..g] do
    sW := tr(p^k, W) - (k ge 2 select p*tr(p^(k-2), W) else 0);
    if IsOdd(k) then
      sW2 := tr(p^k, W2) - (k ge 2 select p*tr(p^(k-2), W2) else 0);
      Append(~S, 2*sW2 - sW);
    else Append(~S, sW); end if;
  end for;
  // S[k] = power sums of Frobenius-of-X_h eigenvalues; build char poly
  c := [-S[1]];
  for j in [2..g] do Append(~c, -(S[j] + &+[c[i]*S[j-i] : i in [1..j-1]]) / j); end for;
  for j in [g+1..2*g-1] do Append(~c, p^(j-g)*c[2*g-j]); end for; Append(~c, p^g);
  key := "[1," cat Join([IntegerToString(Integers()!x) : x in c], ",") cat "]";
  printf "XCHECK %o %o W=%o w%o p=%o ES=%o expected=%o match=%o\n", D, N, f[4], Q, p, key, f[7], key eq f[7];
end for;
quit;
