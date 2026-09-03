// Does the Im(z0) distribution explain the DOCUMENTED base-dependent precision?
//
// VectorValuedForm.m records, at the same Prec := 80:
//     58_5   reldiff <= 5.9e-33   (33 digits)   M = 580
//     34_11  reldiff <= 1.07e-18  (18 digits)   M = 748
//     10_61  reldiff ~ 1          (catastrophic) M = 1220
// If longer coset words drive z0 toward the real axis as M grows, that is the MECHANISM behind
// "achievable precision is base-dependent" -- currently recorded as an observation with no cause.
AttachSpec("ShimuraQuotients.spec");

CC := ComplexField(80); ii := CC.1;
tau0 := CC!0.31 + CC!1.31*ii;
tau1 := CC!(-0.57) + CC!1.73*ii;

slashdata := function(word, tau)
    z := tau; factor := CC!1;
    for i := #word to 1 by -1 do
        if word[i][1] eq "S" then factor /:= Sqrt(z); z := -1/z;
        else z := z + word[i][2]; end if;
    end for;
    return factor, z;
end function;

RR := RealField(6);
printf "%-9o %-6o %-7o %-7o %-13o %-13o %-13o %o\n",
       "base","M","#words","maxlen","min Im(z0)","median","max","known precision";

cases := [ <15,2,"exact (test)">, <58,5,"33 digits">, <34,11,"18 digits">, <10,61,"CATASTROPHIC"> ];
for c in cases do
    D := c[1]; N := c[2];
    M := IsOdd(D*N) select 4*D*N else 2*D*N;
    reps := VVCosetReps(M);
    words := [ VVSTWord(g) : g in reps ];
    im0 := [];
    for w in words do
        _, z0 := slashdata(w, tau0);
        Append(~im0, Imaginary(z0));
    end for;
    s := Sort(im0);
    maxlen := Maximum([ #w : w in words ]);
    printf "%-9o %-6o %-7o %-7o %-13o %-13o %-13o %o\n",
        Sprintf("%o_%o",D,N), M, #words, maxlen,
        RR!s[1], RR!s[#s div 2], RR!s[#s], c[3];
end for;
quit;
