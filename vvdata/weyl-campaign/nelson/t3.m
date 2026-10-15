SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/nelson/nelson.m";

D := 6; N := 1; k := 4; prec := 30;
CC := ComplexField(prec); RR := RealField(prec); I := CC.1;
cut := RR ! 35;

f := Newforms(CuspForms(Gamma0(D*N), k))[1][1];
nmax := 3000;
q := qExpansion(f, nmax+1);
afs := [ Rationals() ! Coefficient(q, n) : n in [1..nmax] ];
af := func< n | n le nmax select CC ! afs[n] else 0 >;

dat := NelsonData(D, N, k, prec);
A := dat`A; O := dat`O; mrep := dat`mrep;

function Act(M, z) return (M[1][1]*z + M[1][2])/(M[2][1]*z + M[2][2]); end function;

z1 := CC ! (0.21 + 1.17*I);
z2 := CC ! (-0.33 + 0.85*I);

t0 := Cputime();
S12 := NelsonPair(dat, z1, z2, af, cut, CC);
S21 := NelsonPair(dat, z2, z1, af, cut, CC);
S11 := NelsonPair(dat, z1, z1, af, cut, CC);
S22 := NelsonPair(dat, z2, z2, af, cut, CC);
printf "time for 4 pair evaluations: %o s\n", Cputime(t0);
printf "S(z1,z2)            = %o\n", S12;
printf "conj S(z2,z1)       = %o\n", ComplexConjugate(S21);
printf "Hermitian residual  = %o\n", Abs(S12 - ComplexConjugate(S21))/Abs(S12);
printf "Cauchy-Schwarz |S12|^2 <= S11*S22 : %o vs %o\n", Abs(S12)^2, S11*S22;

// --- modular transformation law.  f_B(gz) = (cz+d)^k f_B(z) for g in Gamma' = iota(O^1),
//     so S(z1, g z2)/S(z1, z2) must equal (c z2 + d)^k / |c z2 + d|^k, of modulus one.
G := FuchsianGroup(O);
gens := [ A ! Quaternion(g) : g in Generators(G) ];
printf "\n#generators of Gamma': %o, reduced norms %o\n", #gens, [Norm(g) : g in gens];
for g in gens do
    if Norm(g) ne 1 then continue; end if;
    Mg := ChangeRing(mrep(g), CC);
    gz2 := Act(Mg, z2);
    if Imaginary(gz2) le 0 then continue; end if;
    Sg := NelsonPair(dat, z1, gz2, af, cut, CC);
    jj := Mg[2][1]*z2 + Mg[2][2];
    pred := jj^k / Abs(jj)^k;
    got := Sg / S12;
    printf "  g: |pred|=%o  got/pred = %o   rel err %o\n", Abs(pred), got/pred, Abs(got/pred - 1);
end for;
quit;
