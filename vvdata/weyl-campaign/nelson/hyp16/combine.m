// Concatenate two hypeval dumps of the same base, cutoff, precision, rho and circle size (made
// against the same Eichler order and z1) into one, so extra circles can be added to a run.
//   magma -b IN1:=ev_6_185_c46.m IN2:=ev_6_185_c46x.m OUT:=ev_6_185_c46all.m combine.m
SetColumns(0);
function Rd(fn)
    return eval (Read(fn) cat
        "\nreturn <DD, NN, WW, GG, PRECV, CUTV, RHO, NCIRC, NPC, HALLS, ALSIG, PTS, OWNER, ISCTR, VALS>;");
end function;
a := Rd(IN1); b := Rd(IN2);
for i in [1,2,3,4,5,6,7,9,10,11] do
    error if a[i] ne b[i], Sprintf("dumps differ in header field %o", i);
end for;
nc1 := a[8]; nc2 := b[8];
PTS := a[12] cat b[12]; OWNER := a[13] cat [ o + nc1 : o in b[13] ]; ISCTR := a[14] cat b[14];
VALS := a[15] cat b[15];
error if #PTS ne (nc1+nc2)*(a[9]+1) or #VALS ne #PTS, "row counts do not tile";
F := Open(OUT, "w");
fprintf F, "// combined from %o and %o\n", IN1, IN2;
fprintf F, "DD := %o; NN := %o; WW := %o; GG := %o; PRECV := %o; CUTV := %o;\n", a[1], a[2], a[3], a[4], a[5], a[6];
fprintf F, "RHO := %o; NCIRC := %o; NPC := %o;\n", a[7], nc1+nc2, a[9];
fprintf F, "HALLS := %o;\nALSIG := %o;\n", a[10], a[11];
fprintf F, "PTS := [\n";
for j in [1..#PTS] do fprintf F, "  [%o, %o]%o\n", PTS[j][1], PTS[j][2], j lt #PTS select "," else ""; end for;
fprintf F, "];\nOWNER := %o;\nISCTR := %o;\n", OWNER, ISCTR;
fprintf F, "VALS := [\n";
g := a[4];
for j in [1..#VALS] do
    fprintf F, "  [%o]%o\n",
        &cat[ Sprintf("[%o, %o]%o", VALS[j][i][1], VALS[j][i][2], i lt g select ", " else "") : i in [1..g] ],
        j lt #VALS select "," else "";
end for;
fprintf F, "];\n";
delete F;
printf "wrote %o: %o circles, %o points\n", OUT, nc1+nc2, #PTS;
print "COMBINE DONE";
quit;
