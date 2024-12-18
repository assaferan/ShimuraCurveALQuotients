declare type ShimuraQuot;

declare attributes ShimuraQuot: D, N, W, g, CurveID, CoveredBy, Covers, IsP1, IsEC, IsHyp, IsSubhyp;

intrinsic Print(x::ShimuraQuot, L::MonStgElt)
    {Print X at level L}
    if L eq "Default" then
        printf "Shimura quotient of level %o, and discriminant %o by Atkin-Lehners %o", x`N, x`D, x`W;
    elif L eq "Magma" then
        keys := ["D", "N", "W", "g", "CurveID", "CoveredBy", "Covers", "IsP1", "IsEC", "IsHyp", "IsSubhyp"];
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
  for a in ["D", "N", "W", "g", "CurveID", "CoveredBy", "Covers", "IsP1", "IsEC", "IsHyp", "IsSubhyp"] do 
    i := Index(keys,a);
    if i ne 0 then
        x``a := data[i][2];
    end if;
  end for;
  return x;
end intrinsic;

intrinsic 'eq'(crv1::ShimuraQuot, crv2::ShimuraQuot) ->BoolElt
    {returns true if they are equal}
    if (crv1`D ne crv2`D) then
    return false;
    end if;
    if (crv1`N ne crv2`N) then
    return false;
    end if;
    if (crv1`W ne crv2`W) then
    return false;
    end if;
    if (crv1`g ne crv2`g) then
    return false;
    end if;
    if (assigned crv1`IsP1) then
    if not assigned crv2`IsP1 then
        return false;
    end if;
    if (crv1`IsP1 ne crv2`IsP1) then
        return false;
    end if;
    end if;
    if (assigned crv1`IsEC) then
    if not assigned crv2`IsEC then
        return false;
    end if;
    if (crv1`IsEC ne crv2`IsEC) then
        return false;
    end if;
    end if;
    if (assigned crv1`IsHyp) then
    if not assigned crv2`IsHyp then
        return false;
    end if;
    if (crv1`IsHyp ne crv2`IsHyp) then
        return false;
    end if;
    end if;
    if (assigned crv1`IsSubhyp) then
    if not assigned crv2`IsSubhyp then
        return false;
    end if;
    if (crv1`IsSubhyp ne crv2`IsSubhyp) then
        return false;
    end if;
    end if;
    return true;
end intrinsic;

