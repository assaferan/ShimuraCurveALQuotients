// REPAIR track, item 2: what distinguishes wi = 2 and wi = 1221 at 10_61?
//
// Hypothesis: this is NOT a malformed form. In M0MultiplierExact the slash constant is
// evaluated at TWO HARDCODED points (VectorValuedForm.m:419-420)
//     tau0 = 0.31 + 1.31i,  tau1 = -0.57 + 1.73i
// and slashdata(w, tau) walks the word applying z -> -1/z and z -> z+n. A word can drive z
// close to the real axis, where DedekindEta(z) blows up -- so k0 = num0/sfun(tau0) can
// overflow while k1 stays sane. That would be an EVALUATION-POINT defect, fixable by choosing
// tau per word, not an arithmetic defect in the base.
//
// words depends ONLY on M = 2*D*N (VectorValuedForm.m:393-394), so this needs no Borcherds
// forms and no pipeline: 10_61 -> M = 1220.
AttachSpec("ShimuraQuotients.spec");

D := 10; N := 61;
M := IsOdd(D*N) select 4*D*N else 2*D*N;
printf "D=%o N=%o  ->  M = %o\n", D, N, M;

CC := ComplexField(80); ii := CC.1;
tau0 := CC!0.31 + CC!1.31*ii;
tau1 := CC!(-0.57) + CC!1.73*ii;

// verbatim from VectorValuedForm.m:410-417
slashdata := function(word, tau)
    z := tau; factor := CC!1;
    for i := #word to 1 by -1 do
        if word[i][1] eq "S" then factor /:= Sqrt(z); z := -1/z;
        else z := z + word[i][2]; end if;
    end for;
    return factor, z;
end function;

t := Realtime();
reps := VVCosetReps(M);
words := [ VVSTWord(g) : g in reps ];
printf "#words = %o  (%os)\n\n", #words, Realtime()-t;

im0 := []; im1 := [];
for wi in [1..#words] do
    _, z0 := slashdata(words[wi], tau0);
    _, z1 := slashdata(words[wi], tau1);
    Append(~im0, Imaginary(z0));
    Append(~im1, Imaginary(z1));
end for;

RR := RealField(6);
printf "%-8o %-14o %-14o %o\n", "wi", "Im(z0)", "Im(z1)", "note";
for wi in [2, 1221] do
    if wi le #words then
        printf "%-8o %-14o %-14o %o\n", wi, RR!im0[wi], RR!im1[wi], "<<< THE TWO FAILING COSETS";
    end if;
end for;
printf "\n";

// Where do the two sit in the overall distribution?
s0 := Sort(im0); s1 := Sort(im1);
printf "Im(z0) over all words: min %o  median %o  max %o\n", RR!s0[1], RR!s0[#s0 div 2], RR!s0[#s0];
printf "Im(z1) over all words: min %o  median %o  max %o\n", RR!s1[1], RR!s1[#s1 div 2], RR!s1[#s1];
printf "\n";

// Rank the words by how small Im(z0) is -- if the hypothesis holds, 2 and 1221 are extremal.
idx0 := [1..#words];
Sort(~idx0, func< a, b | im0[a] lt im0[b] select -1 else (im0[a] gt im0[b] select 1 else 0) >);
printf "10 SMALLEST Im(z0) (wi : Im(z0) / Im(z1)):\n";
for j in [1..Minimum(10, #idx0)] do
    wi := idx0[j];
    printf "   wi %-6o Im(z0) %-14o Im(z1) %o%o\n", wi, RR!im0[wi], RR!im1[wi],
        (wi in [2,1221]) select "   <<<" else "";
end for;
printf "\n";
idx1 := [1..#words];
Sort(~idx1, func< a, b | im1[a] lt im1[b] select -1 else (im1[a] gt im1[b] select 1 else 0) >);
printf "10 SMALLEST Im(z1) (the control -- k1 was SANE at both failing cosets):\n";
for j in [1..Minimum(10, #idx1)] do
    wi := idx1[j];
    printf "   wi %-6o Im(z1) %-14o Im(z0) %o%o\n", wi, RR!im1[wi], RR!im0[wi],
        (wi in [2,1221]) select "   <<<" else "";
end for;

// How bad is an eta at these heights?  eta(z) ~ q^(1/24), q = e(z); small Im -> huge under 1/eta.
printf "\nrank of the two failing cosets by ascending Im(z0): ";
for wi in [2,1221] do
    p := Index(idx0, wi);
    printf "wi %o -> rank %o of %o;  ", wi, p, #words;
end for;
printf "\n";
quit;
