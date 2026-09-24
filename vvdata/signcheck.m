// Reproduce the X0^14(5) failure against the fixed routine, on the exact values from the crash.
AttachSpec("ShimuraQuotients.spec");
s      := [* Infinity(), 5/8, 55/128, 0, 1, 1/4, 1/8 *];
stilde := [* Infinity(), 35/8, 585/128, 1, 0, 3/4, 9/4 *];
ds     := [ -4, -35, -280, -11, -91, -84, -51 ];
degs   := [ 1, 1, 1, 1, 1, 1, 2 ];
scale_tilde := stilde[Index(s,0)];  scale := s[Index(stilde,0)];
printf "scale = %o, scale_tilde = %o\n", scale, scale_tilde;
inf_zero := [Index(s,0), Index(stilde,0), Index(s,Infinity())];
for i in [1..#s] do
    if i in inf_zero or degs[i] ne 1 then
        printf "  d = %-6o SKIPPED (%o)\n", ds[i],
               i in inf_zero select "normalising point" else "degree 2";
        continue;
    end if;
    sols := [ [e1,e2] : e1,e2 in [-1,1] | e1*s[i]/scale + e2*stilde[i]/scale_tilde eq 1 ];
    N := 5; df := FundamentalDiscriminant(ds[i]);
    printf "  d = %-6o %-12o sign choices = %o  %o\n", ds[i],
           (df mod N eq 0) select "NON-firing" else "FIRING", sols,
           IsEmpty(sols) select "  <== INCONSISTENT" else "";
end for;
quit;
