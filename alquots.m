declare type ShimuraQuot;

declare attributes ShimuraQuot: D, N, W, g, CurveID, CoveredBy, Covers, IsP1, IsEC, IsHyp, IsSubhyp;

intrinsic Print(x::ShimuraQuot, L::MonStgElt)
    {Print X at level L}
    if L eq "Default" then
        printf "Shimura quotient of level %o, and discriminant %o by Atkin-Lehners %o", x`N, x`D, x`W;
    elif L eq "Magma" then
        keys := ["D", "N", "W", "g", "CurveID", "CoveredBy", "Covers", "IsP1", "IsEC", "IsHyp", "IsSubhyp"];
        data := [* <k, x``k> : k in keys | assigned x``k *];
        printf "CreateShimuraQuot(%m);\n", data;
    end if;

end intrinsic;

intrinsic CreateShimuraQuot(D ::RngIntElt,N ::RngIntElt,W ::SetEnum) -> ShimuraQuot
  {Create Shimura Quotient with disc D, level N, and ALs W}
  x := New(ShimuraQuot);
  assert GCD(D,N) eq 1;
  x`D := D;
  x`N := N;
  x`W := W;
  return x;
end intrinsic;

intrinsic CreateShimuraQuot(data :: List) -> ShimuraQuot
  {Create Shimura Quotient from the magma printout}
  x := New(ShimuraQuot);
  keys := [d[1] : d in data];
  for a in ["D", "N", "W", "g", "CurveID", "CoveredBy", "Covers", "IsP1", "IsEC", "IsHyp", "IsSubhyp"] do 
    i := Index(keys,a);
    if i ne 0 then
        x``a := data[i][2];
    end if;
  end for;
  return x;
end intrinsic;


