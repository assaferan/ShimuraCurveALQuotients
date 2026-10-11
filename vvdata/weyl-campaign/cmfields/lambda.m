// Lambda = O cap {x : x iota = -iota x} for optimal embeddings iota of the order R of disc d into O.
// Claim under test: -n|_Lambda / m is a primitive form of discriminant d whose class in
// Pic(R)/Pic(R)^2 is the conjugation class [a] of GR Lemma 5.10 (m = D_R N*_R).
function anticommuting_form(O, iota, m)
    B := Algebra(O); bas := Basis(O);
    M := Matrix(Rationals(), 4, 4, [Eltseq(B!(iota*bas[a] + bas[a]*iota)) : a in [1..4]]);
    Kr := Kernel(M);                                   // coordinates w.r.t. bas of V = anticommuting plane
    rows := [];
    for v in Basis(Kr) do
        den := LCM([Denominator(c) : c in Eltseq(v)]);
        Append(~rows, [Integers()| c*den : c in Eltseq(v)]);
    end for;
    S := Saturation(Matrix(Integers(), rows));         // Z-basis of O cap V
    Vs := [&+[S[r][k]*bas[k] : k in [1..4]] : r in [1..2]];
    G := Matrix(Rationals(), 2, 2, [Trace(Vs[a]*Conjugate(Vs[b]))/2 : a, b in [1..2]]);
    return -G[1][1]/m, -2*G[1][2]/m, -G[2][2]/m;
end function;
procedure report(D, N, d, m)
    B := QuaternionAlgebra(D);
    O := N eq 1 select MaximalOrder(B) else Order(MaximalOrder(B), N);
    bas := Basis(O);
    // generator of R: sqrt(d) if d = 0 mod 4 (x^2 = d, take x = sqrt(d)), else (1+sqrt(d))/2 has trace 1 norm (1-d)/4
    seen := {}; nemb := 0;
    for c1, c2, c3, c4 in [-10..10] do
        x := c1*bas[1] + c2*bas[2] + c3*bas[3] + c4*bas[4];
        if d mod 4 eq 0 then
            if Trace(x) ne 0 or Norm(x) ne -d div 4 then continue; end if;   // x = sqrt(d/4), the generator of the order of disc d
            // optimal: Z[x] = O cap Q(x): no (u + x)/p in O for p | conductor-part
            ok := true;
            for p in PrimeDivisors(Integers()!(-d)) do
                for u in [0..p-1] do if (u + x)/p in O then ok := false; end if; end for;
            end for;
            if not ok then continue; end if;
            iota := x;
        else
            if Trace(x) ne 1 or Norm(x) ne (1 - d) div 4 then continue; end if;   // x = (1+sqrt d)/2
            ok := true;
            for p in PrimeDivisors(Integers()!(-d)) do
                for u in [0..p-1] do if (u + x)/p in O then ok := false; end if; end for;
            end for;
            if not ok then continue; end if;
            iota := 2*x - 1;
        end if;
        nemb +:= 1;
        a, b, c := anticommuting_form(O, iota, m);
        if not (IsIntegral(a) and IsIntegral(b) and IsIntegral(c)) then printf "  iota=%o: form (%o,%o,%o) NOT integral\n", iota, a, b, c; continue; end if;
        disc := Integers()!(b^2 - 4*a*c);
        key := <disc, ReducedForm(QuadraticForms(disc)![Integers()|a, b, c])>;
        if key in seen then continue; end if; Include(~seen, key);
        printf "  -n|_Lambda/%o = (%o, %o, %o): disc %o, reduced %o, content %o\n", m, a, b, c, disc, key[2], GCD([Integers()|a,b,c]);
    end for;
    printf "  (%o embeddings found in the box)\n", nemb;
end procedure;
D := StringToInteger(D); N := StringToInteger(N); d := StringToInteger(d); m := StringToInteger(m);
printf "D=%o, N=%o, d=%o, m=%o:\n", D, N, d, m;
report(D, N, d, m);
