// Noether's conic on subsets of the rows of a dump: per chunk of rows and overall.
//   magma -b IN:=ev_55_7_c40.m CH:=51 conicdiag.m < /dev/null
SetColumns(0);
dmp := eval (Read(IN) cat
   "\nreturn <DD, NN, WW, GG, PRECV, CUTV, RHO, NCIRC, NPC, HALLS, ALSIG, PTS, OWNER, ISCTR, VALS>;");
PRECV := dmp[5]; PTS := dmp[12]; VALS := dmp[15]; g := dmp[4]; ISCTR := dmp[14]; NPC := dmp[9];
CC := ComplexField(PRECV); RR := RealField(PRECV); Ic := CC.1; R6 := RealField(6);
pts := [ CC | PTS[j][1] + Ic*PTS[j][2] : j in [1..#PTS] ];
Gv := [ [ CC | (VALS[j][i][1] + Ic*VALS[j][i][2]) / Imaginary(pts[j]) : i in [1..g] ] : j in [1..#PTS] ];
mons := [ [a,b] : b in [a..g], a in [1..g] ];
function ConicSV(idx)
    Mx := Matrix(CC, #idx, #mons, &cat[ [ Gv[j][m[1]]*Gv[j][m[2]] : m in mons ] : j in idx ]);
    for i in [1..#idx] do
        nm := Sqrt(&+[ Abs(Mx[i][c])^2 : c in [1..#mons] ]);
        for c in [1..#mons] do Mx[i][c] := Mx[i][c]/nm; end for;
    end for;
    n2 := #idx;
    MR := Matrix(RR, 2*n2, 2*#mons, [ RR ! 0 : i in [1..4*n2*#mons] ]);
    for i in [1..n2] do for j in [1..#mons] do
        a := Real(Mx[i][j]); b := Imaginary(Mx[i][j]);
        MR[i][j] := a; MR[i][#mons+j] := -b; MR[n2+i][j] := b; MR[n2+i][#mons+j] := a;
    end for; end for;
    Sg := SingularValueDecomposition(MR);
    rr := Min(Nrows(MR), Ncols(MR));
    sv := Reverse([ Abs(Sg[i][i]) : i in [1..rr] ]);
    return [ R6 ! sv[i] : i in [1..#sv by 2] ];
end function;
ch := StringToInteger(CH);
printf "all centres: %o\n", ConicSV([ j : j in [1..#PTS] | ISCTR[j] eq 1 ]);
printf "all rows with j mod %o eq 3: %o\n", NPC, ConicSV([ j : j in [1..#PTS] | (j mod NPC) eq 3 ]);
s := 1;
while s le #PTS do
    e := Min(s + ch - 1, #PTS);
    idx := [ j : j in [s..e] | (j mod 4) eq 1 ];
    printf "rows %o..%o (every 4th): %o\n", s, e, ConicSV(idx);
    s := e + 1;
end while;
// swapping columns 2 and 3 on one chunk at a time
for k in [1..Ceiling(#PTS/ch)] do
    G2 := Gv;
    for j in [(k-1)*ch+1 .. Min(k*ch, #PTS)] do G2[j] := [ Gv[j][1], Gv[j][3], Gv[j][2] ]; end for;
    Gv := G2;
    printf "columns 2,3 swapped on chunk %o: centres %o\n", k, ConicSV([ j : j in [1..#PTS] | ISCTR[j] eq 1 ]);
    Gv := dmp[15]; Gv := [ [ CC | (VALS[j][i][1] + Ic*VALS[j][i][2]) / Imaginary(pts[j]) : i in [1..g] ] : j in [1..#PTS] ];
end for;
print "CONICDIAG DONE";
quit;
