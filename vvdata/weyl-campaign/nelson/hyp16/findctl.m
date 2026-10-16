// Find more negative controls IN THE STRUCTURAL CLASS OF THE SIXTEEN: genus-3 quotients with NO
// Atkin-Lehner overgroup of genus 0, so that a hyperelliptic involution there would have to be
// exceptional, and whose hyperellipticity is ALREADY KNOWN exactly because the repository stores
// a function field for them.  Those are the curves on which a false positive would show.
// Also requires that the newforms of level exactly D*N span, or Theorem 3.1 cannot be used.
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
dir := "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/data/models/";
fs := Split(Pipe("ls " cat dir, ""), "\n");
out := [];
for fn in fs do
    if #fn lt 10 or fn[1..7] ne "models_" then continue; end if;
    bits := Split(fn[8..#fn-2], "_");
    if #bits ne 2 then continue; end if;
    D := StringToInteger(bits[1]); N := StringToInteger(bits[2]);
    if not (IsSquarefree(D*N) and GCD(D,N) eq 1) then continue; end if;
    if D*N gt 700 then continue; end if;             // keep the evaluation cheap
    models := eval (Read(dir cat fn) cat "\nreturn models;");
    halls := [ m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1 ];
    for key in Keys(models) do
        W := Set(key);
        g := GenusShimuraCurveQuotient(D, N, W);
        if g ne 3 then continue; end if;
        // no Atkin-Lehner overgroup of genus 0  ->  any hyperelliptic involution is exceptional
        bad := false;
        for m in halls do
            if m in W then continue; end if;
            WW := { (a*m) div GCD(a,m)^2 : a in W } join W;
            if #WW ne 2*#W then continue; end if;
            if GenusShimuraCurveQuotient(D, N, WW) eq 0 then bad := true; break; end if;
        end for;
        if bad then continue; end if;
        if #models[key] eq 0 then continue; end if;
        Append(~out, <D, N, Sort(SetToSequence(W)), D*N, models[key][1][2]>);
    end for;
end for;
printf "genus-3 stored quotients with NO genus-0 Atkin-Lehner overgroup, D*N <= 700:\n";
Sort(~out, func< a,b | a[4]-b[4] >);
for t in out do
    D := t[1]; N := t[2]; W := Set(t[3]);
    // does Theorem 3.1 reach the whole space?
    nmax := 800;
    NFs := Newforms(CuspForms(Gamma0(D*N), 2));
    cnt := 0;
    for orb in NFs do for f in orb do
        if Level(Parent(f)) ne D*N then continue; end if;
        KF := BaseRing(Parent(f)); qq := qExpansion(f, nmax+1);
        cj := Degree(KF) eq 1 select [ ComplexField(30) | 0 ]
              else Conjugates(KF.1 : Precision := 30);
        idx := Degree(KF) eq 1 select 1 else Index(orb, f);
        cs := [ Degree(KF) eq 1 select ComplexField(30)!(Rationals()!Coefficient(qq,n))
                else Evaluate(Polynomial(Rationals(), Eltseq(Coefficient(qq,n))), cj[idx])
                : n in [1..nmax] ];
        inv := true;
        for m in W do
            if m eq 1 then continue; end if;
            e := ComplexField(30) ! 1;
            for p in PrimeDivisors(m) do
                e *:= (D mod p eq 0) select cs[p] else -cs[p];
            end for;
            if Abs(e - 1) gt 0.01 then inv := false; end if;
        end for;
        if inv then cnt +:= 1; end if;
    end for; end for;
    printf "  %o_%o  W = %o  D*N = %o  newforms of exact level %o/3 %o  shape %o\n",
       D, N, t[3], t[4], cnt, cnt eq 3 select "<= USABLE" else "(short)",
       Type(t[5]) eq MonStgElt select "CRV" else "poly";
end for;
printf "FINDCTL DONE\n";
quit;
