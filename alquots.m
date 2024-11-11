declare type ShimuraQuot;

declare attributes ShimuraQuot: D, N, W, g, CurveID, CoveredBy, Covers, IsP1, IsEC, IsHyp, IsSubhyp;

intrinsic Print(x::ShimuraQuot)
  {Print x}
  printf "Shimura quotient of level %o, and discriminant %o by Atkin-Lehners %o", x`N, x`D, x`W;
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