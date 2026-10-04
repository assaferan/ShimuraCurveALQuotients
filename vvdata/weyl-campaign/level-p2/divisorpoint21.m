// A second gap for the -log m rule at a point on the divisor: X_0^21(2), d = -4 (s = oo), where the
// forms' poles come from x = lambda_0 (m = 1), x = 3 lambda_0 (m = 9) and the cusp-0 pair (m = 1/4).
// Truth: |f_k| = C_k prod |s - s_i|^{m_i} with the Guo-Yang s-values; at s = oo a quotient of forms
// with equal total degree tends to the ratio of the C_k.  C_k is read at d = -7 (s = -7).
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
D := 21; N := 2;
sv := AssociativeArray();
for t in [<-7,-7>,<-15,-5/3>,<-16,1>,<-28,1/9>,<-60,9>,<-84,-3>,<-100,1/5>,<-112,25>,<-120,-1/3>,<-148,37/9>,<-168,0>,<-228,-25/3>,<-232,-32>,<-280,-35/9>,<-312,-8/3>,<-372,-3/4>,<-408,-75>,<-420,21>,<-532,-19/4>,<-708,25/48>,<-840,-16/3>] do sv[t[1]] := Rationals()!t[2]; end for;
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
ks := Sort([k : k in Keys(fs)]);
Ld := ShimuraCurveLattice(D, N);
function cell_value(x)
    v := Rationals()!1;
    for q in Keys(x`log_coeffs) do c := x`log_coeffs[q]; if not IsIntegral(c) then return false; end if; v *:= (Rationals()!q)^(Integers()!c); end for;
    return v;
end function;
// base points for C_k: the first table point off the form's divisor and coprime to the level
base := [-7, -15, -16, -28, -60, -84, -100, -112, -120, -148, -168, -228, -232, -280, -312, -372, -408];
V := AssociativeArray();
for d in Keys(sv) do if GCD(d, N) eq 1 then V[d] := SchoferFormula([fs[k] : k in ks], d, D, N, Ld); end if; end for;
valsat := func< d | V[d] >;
// forms 12..15 have half-integral log coefficients (their principal parts are not integral, so they
// are outside Theorem B's hypothesis); their DOUBLES are inputs, and Schofer's formula is linear in F.
// Work throughout with 2 f: C2[k] = C_k^2, degrees as for f.
C2 := AssociativeArray(); degs := AssociativeArray(); dvs := AssociativeArray();
for i->k in ks do
    dv := DivisorOfBorcherdsForm(fs[k], star); dvs[k] := dv;
    for pr in dv do error if pr[1] ne -4 and not IsDefined(sv, pr[1]), Sprintf("form %o: divisor point %o not in the table", k, pr[1]); end for;
    model := func< d | &*[Rationals() | AbsoluteValue(sv[d] - sv[pr[1]])^(Integers()!(2*pr[2])) : pr in dv | pr[1] ne -4] >;
    ok := false;
    for d0 in base do
        if GCD(d0, N) ne 1 or exists{pr : pr in dv | pr[1] eq d0} then continue; end if;
        cv := cell_value(2*valsat(d0)[i]);
        if cv cmpeq false then printf "form %o at %o: double irrational %o\n", k, d0, valsat(d0)[i]; continue; end if;
        C2[k] := cv / model(d0); ok := true; printf "form %-3o C^2 read at %o: %o\n", k, d0, C2[k]; break;
    end for;
    error if not ok, Sprintf("form %o: no base point", k);
    degs[k] := &+[Rationals() | pr[2] : pr in dv | pr[1] ne -4];
    foo := qExpansionAtoo(fs[k], 1); f0 := qExpansionAt0(fs[k], 1);
    printf "form %-3o divisor %o  deg %o  poles at oo %o, at 0 (q^(1/168)) %o\n", k, dv, degs[k],
        [m : m in [1..-Valuation(foo)] | Coefficient(foo,-m) ne 0], [m : m in [1..-Valuation(f0)] | Coefficient(f0,-m) ne 0];
end for;
nch := 0;
for d in Keys(sv) do
    if GCD(d, N) ne 1 then continue; end if;       // d = -420 is not coprime to the level: the known p | gcd(d, N) gap
    vals := valsat(d);
    for i->k in ks do
        dv := dvs[k];
        if exists{pr : pr in dv | pr[1] eq d} then continue; end if;
        model := &*[Rationals() | AbsoluteValue(sv[d] - sv[pr[1]])^(Integers()!(2*pr[2])) : pr in dv | pr[1] ne -4];
        cv := cell_value(2*vals[i]);
        if cv cmpeq false then printf "form %o at d = %o: double irrational: %o\n", k, d, vals[i]; continue; end if;
        if cv ne C2[k]*model then printf "MISMATCH form %o at d = %o: %o vs %o\n", k, d, cv, C2[k]*model; else nch +:= 1; end if;
    end for;
end for;
printf "%o off-divisor values (of 2f) agree with the divisor model\n", nch;
SetVerbose("ShimuraQuotients", 1);
v4 := SchoferFormula([fs[k] : k in ks], -4, D, N, Ld);
SetVerbose("ShimuraQuotients", 0);
printf "\nat d = -4 (s = oo): (2 f_a)^deg_b / (2 f_b)^deg_a  must equal (C_a^2)^deg_b / (C_b^2)^deg_a\n";
for i->a in ks, j->b in ks do
    if a ge b or degs[a] eq 0 or degs[b] eq 0 then continue; end if;
    na := Integers()!(2*degs[b]); nb := Integers()!(2*degs[a]);
    got := cell_value(na*v4[i] - nb*v4[j]); want := C2[a]^(na div 2) / C2[b]^(nb div 2);
    if got cmpeq false then printf "  (%o,%o) irrational: %o\n", a, b, na*v4[i] - nb*v4[j]; continue; end if;
    printf "  f_%o^%o / f_%o^%o = %o   divisors give %o   %o\n", a, na, b, nb, got, want, got eq want select "ok" else "MISMATCH";
end for;
