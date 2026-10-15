// The other half of the cost, and the one that decides the reach: where do the eigenvalues come
// from?  Theorem 3.1 needs a_f(n) for every n up to DN*(cut/4pi)^2, and the lattice enumeration
// is useless without them.  This times the classical route -- decompose S_2(Gamma_0(DN)) into
// newforms and expand one of them -- at increasing level, so that the point at which it stops
// being the right tool is measured rather than guessed.
//
//   magma -b LEVELS:=282,385,429,510 CUT:=32 eig.m < /dev/null
SetColumns(0);
levs := [ StringToInteger(x) : x in Split(LEVELS, ",") ];
cutv := StringToInteger(CUT);
RR := RealField(20);
C2 := (RR ! cutv / (4*Pi(RR)))^2;
printf "%-8o %-10o %-8o %-12o %-12o %-12o %o\n",
   "level", "dim S_2", "n_max", "newforms s", "expand s", "orbits", "degrees";
for L in levs do
    nmax := Ceiling(L*C2) + 10;
    t0 := Cputime();
    S := CuspForms(Gamma0(L), 2);
    d := Dimension(S);
    NF := Newforms(S);
    tnf := Cputime(t0);
    t1 := Cputime();
    f := NF[#NF][1];
    _ := qExpansion(f, nmax+1);
    tq := Cputime(t1);
    printf "%-8o %-10o %-8o %-12o %-12o %-12o %o\n", L, d, nmax,
        RealField(5)!tnf, RealField(5)!tq, #NF,
        [ Degree(BaseRing(Parent(o[1]))) : o in NF ];
end for;
printf "EIG DONE\n";
quit;
