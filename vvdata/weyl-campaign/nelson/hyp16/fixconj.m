// Repair a chunked dump whose chunks disagree on which Galois conjugate of a degree-2 newform
// orbit is which column.  Magma's generator of the coefficient field, hence the order of its
// conjugates, is not the same from one process to the next, so chunks computed separately can
// have the two conjugate columns swapped relative to each other; each chunk alone then lies on
// a canonical conic and the merged set on none.  For every pattern of swaps over the chunks
// (chunk 1 fixed) the conic residual of the centres is measured, and the pattern with the
// smallest residual is written out if it is clean.
//   magma -b IN:=ev_6_323_c46.m OUT:=ev_6_323_c46f.m CH:=59 C1:=2 C2:=3 fixconj.m < /dev/null
SetColumns(0);
dmp := eval (Read(IN) cat
   "\nreturn <DD, NN, WW, GG, PRECV, CUTV, RHO, NCIRC, NPC, HALLS, ALSIG, PTS, OWNER, ISCTR, VALS>;");
g := dmp[4]; PRECV := dmp[5]; PTS := dmp[12]; ISCTR := dmp[14]; VALS := dmp[15]; NPC := dmp[9];
CC := ComplexField(PRECV); RR := RealField(PRECV); Ic := CC.1; R6 := RealField(6);
ch := StringToInteger(CH); c1 := StringToInteger(C1); c2 := StringToInteger(C2);
pts := [ CC | PTS[j][1] + Ic*PTS[j][2] : j in [1..#PTS] ];
nchunks := Ceiling(#PTS/ch);
mons := [ [a,b] : b in [a..g], a in [1..g] ];
function ConicRes(V)
    idx := [ j : j in [1..#PTS] | ISCTR[j] eq 1 ];
    Gv := [ [ CC | (V[j][i][1] + Ic*V[j][i][2]) / Imaginary(pts[j]) : i in [1..g] ] : j in idx ];
    Mx := Matrix(CC, #idx, #mons, &cat[ [ Gv[r][m[1]]*Gv[r][m[2]] : m in mons ] : r in [1..#idx] ]);
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
    return sv[1], sv[3]/sv[1];
end function;
function Swapped(V, pattern)
    W := V;
    for k in [2..nchunks] do
        if pattern[k-1] eq 0 then continue; end if;
        for j in [(k-1)*ch+1 .. Min(k*ch, #PTS)] do
            row := V[j]; t := row[c1]; row[c1] := row[c2]; row[c2] := t; W[j] := row;
        end for;
    end for;
    return W;
end function;
best := 0; bestres := RR ! 1; bestgap := RR ! 0;
for bits in CartesianPower([0, 1], nchunks - 1) do
    pat := [ bits[i] : i in [1..nchunks-1] ];
    res, gap := ConicRes(Swapped(VALS, pat));
    printf "swap pattern %o: conic residual %o, gap %o\n", pat, R6 ! res, R6 ! gap;
    if res lt bestres then bestres := res; bestgap := gap; best := pat; end if;
end for;
printf "best pattern %o, residual %o, gap %o\n", best, R6 ! bestres, R6 ! bestgap;
error if bestgap lt 10^6, "no swap pattern gives a clean conic; the inconsistency is something else";
V := Swapped(VALS, best);
txt := Read(IN);
p := Position(txt, "VALS := [");
F := Open(OUT, "w");
fprintf F, "%o", txt[1..p-1];
fprintf F, "// conjugate columns %o, %o swapped on chunks with pattern %o (fixconj.m)\n", c1, c2, best;
fprintf F, "VALS := [\n";
for j in [1..#V] do
    fprintf F, "  [%o]%o\n",
        &cat[ Sprintf("[%o, %o]%o", V[j][i][1], V[j][i][2], i lt g select ", " else "") : i in [1..g] ],
        j lt #V select "," else "";
end for;
fprintf F, "];\n";
delete F;
printf "wrote %o\nFIXCONJ DONE\n", OUT;
quit;
