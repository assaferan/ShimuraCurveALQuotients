// Has Morita's model already been blown down?
//
// At p | D the Weierstrass points with singular specialisation are exactly those by an order in
// Q(sqrt(-n)) with p | n (Thm 4.4), and every singular point of the fibre carries exactly two of
// them (Thm 2.3).  So twice the number of singular points must equal that count.  When it does
// not, the quotient of Morita's model is not yet a Weierstrass model and exceptional components
// have to be blown down, which changes the thicknesses by (4.15).
//
// The implementation does not do the blow-down.  It gets away with it whenever the blow-down
// merges singular points along a single chain, because (4.15) then leaves the SUM of the
// thicknesses -- all that disc R needs -- unchanged.  This script says where to be careful.
//
//   magma -b logfile:=out.txt vvdata/weyl-campaign/molina/diag.m < /dev/null
// Run from the repository root.

AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
load "vvdata/weyl-campaign/molina/molina.m";
L := Open(logfile, "w");
procedure Say(s) fprintf L, "%o\n", s; Flush(L); end procedure;

cases := [<35,5,2>, <51,17,2>, <57,3,2>, <65,13,2>, <65,5,2>, <69,23,2>, <85,5,2>, <85,85,2>,
          <39,1,3>, <55,1,3>, <115,115,2>, <161,161,2>, <115,5,4>];

for c in cases do
    D := c[1]; m := c[2]; g := c[3];
    W := (m eq 1) select {Integers() | 1} else {Integers() | 1, m};
    wd := MolWeierstrassData(D, m, g);
    rows := wd[1][3];
    for p in PrimeDivisors(D) do
        th := MolQuotientThicknesses(D, p, W);
        nsing := &+[Integers() | r[2]*r[3] : r in rows | r[1] mod p eq 0];
        Say(Sprintf("D = %o, w_%o, p = %o: %o Weierstrass points specialise to a singular point,"
            * " the fibre has %o singular points (twice that is %o) -- %o",
            D, m, p, nsing, #th, 2*#th,
            nsing eq 2*#th select "no blow-down" else "BLOWN DOWN"));
    end for;
end for;
Say("DIAG DONE");
quit;
