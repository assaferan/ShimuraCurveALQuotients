// Criterion H: at an odd prime p || N with p !| D, the special fibre of X = X_0^D(N)/W has the
// component Z = X_0^D(N/p)/W', W' = {w in W : p !| w} (Deligne-Rapoport for D = 1; Buzzard and
// Helm for D > 1).  If g(Z) >= 2 and Z/w_u has genus 0 for an Atkin-Lehner u !in W', the
// hyperelliptic involution of X, extended to the stable model, must reduce to w_u, so iota = w_u.

intrinsic ComponentInvolutionNotHyperelliptic(D::RngIntElt, N::RngIntElt, W::SetEnum) -> BoolElt, MonStgElt
{True, with a witness, if some odd p || N with p !| D has a component Z = X_0^D(N/p)/W' of genus
at least 2 and an Atkin-Lehner u, p !| u, u !in W', with Z/w_u of genus 0 and X/w_u of positive
genus; then X = X_0^D(N)/W is not geometrically hyperelliptic.  Uses genera only.}
    star := func<a, b | a*b div GCD(a, b)^2>;
    for p in PrimeDivisors(N) do
        if p eq 2 or Valuation(N, p) ne 1 or D mod p eq 0 then continue; end if;
        M := N div p;
        Wc := {Integers() | w : w in W | w mod p ne 0};
        h := GenusShimuraCurveQuotient(D, M, Wc);
        if h lt 2 then continue; end if;
        for u in [d : d in Divisors(D*M) | GCD(d, (D*M) div d) eq 1 and d notin Wc] do
            if GenusShimuraCurveQuotient(D, M, Wc join {star(u, w) : w in Wc}) ne 0 then continue; end if;
            gX := GenusShimuraCurveQuotient(D, N, W join {star(u, w) : w in W});
            if gX gt 0 then
                return true, Sprintf("ComponentInvolution p=%o case=%o component=X_0^%o(%o)/%o g=%o u=%o g(X/w_u)=%o",
                                     p, Wc eq W select "A" else "B", D, M, Wc, h, u, gX);
            end if;
        end for;
    end for;
    return false, _;
end intrinsic;

intrinsic FilterByComponentInvolution(~curves::SeqEnum)
{Mark the undecided curves of genus at least 3 that ComponentInvolutionNotHyperelliptic rules out.}
    for i->X in curves do
        if assigned X`IsSubhyp then continue; end if;
        if X`g lt 3 then continue; end if;
        ok, witness := ComponentInvolutionNotHyperelliptic(X`D, X`N, X`W);
        if ok then
            curves[i]`IsSubhyp := false;
            curves[i]`IsHyp := false;
            curves[i]`TestInWhichProved := witness;
        end if;
        if (i mod 200 eq 0) then
            vprintf ShimuraQuotients, 1: "i = %o/%o\n", i, #curves;
        end if;
    end for;
end intrinsic;
