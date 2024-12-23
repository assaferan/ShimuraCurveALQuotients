declare type ShimuraQuot;

declare attributes ShimuraQuot: D, N, W, g, CurveID, CoveredBy, Covers, IsP1, IsEC, IsHyp, IsSubhyp, TestInWhichProved;

intrinsic Print(x::ShimuraQuot, L::MonStgElt)
    {Print X at level L}
    if L eq "Default" then
        printf "Shimura quotient of level %o, and discriminant %o by Atkin-Lehners %o", x`N, x`D, x`W;
    elif L eq "Magma" then
        keys := ["D", "N", "W", "g", "CurveID", "CoveredBy", "Covers", "IsP1", "IsEC", "IsHyp", "IsSubhyp", "TestInWhichProved"];
        data := [* <k, x``k> : k in keys | assigned x``k *];
        printf "CreateShimuraQuot(%m);", data;
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
  for a in ["D", "N", "W", "g", "CurveID", "CoveredBy", "Covers", "IsP1", "IsEC", "IsHyp", "IsSubhyp", "TestInWhichProved"] do 
    i := Index(keys,a);
    if i ne 0 then
        x``a := data[i][2];
    end if;
  end for;
  return x;
end intrinsic;

intrinsic 'eq'(crv1::ShimuraQuot, crv2::ShimuraQuot) ->BoolElt
    {returns true if they are equal}

    for attr in ["D", "N", "W", "g", "CurveID"] do
        if (crv1``attr ne crv2``attr) then
            return false;
        end if;
    end for;

    for attr in ["CoveredBy", "Covers", "IsP1", "IsEC", "IsHyp", "IsSubhyp", "TestInWhichProved"] do
        if (assigned crv1``attr) then
            if not assigned crv2``attr then
                return false;
            end if;
            if (crv1``attr ne crv2``attr) then
                return false;
            end if;
        end if;
    end for;
    
    return true;
end intrinsic;