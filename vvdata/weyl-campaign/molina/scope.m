// How far does Molina's method reach among the 43 missing models?  For each one, run Step 1 only
// (cheap: no ring class field is built) and report whether the hypotheses hold and what the
// degrees of the number fields would be -- IndexFormEquation, which Step 4 needs, is implemented
// for degree 3 and 4 only.
//   magma -b logfile:=... scope.m < /dev/null

AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
load "vvdata/weyl-campaign/molina/molina.m";
L := Open(logfile, "w");
procedure Say(s) fprintf L, "%o\n", s; Flush(L); end procedure;

targets := [
  <85,[1,5],2>, <85,[1,85],2>, <91,[1,91],2>, <115,[1,115],2>, <123,[1,123],2>,
  <133,[1,7,19,133],2>, <141,[1,141],2>, <145,[1,5,29,145],2>, <155,[1,155],2>, <161,[1,161],2>,
  <177,[1,3,59,177],2>, <187,[1,11,17,187],2>, <213,[1,3,71,213],2>, <237,[1,79],2>,
  <249,[1,3,83,249],2>, <265,[1,5,53,265],2>, <287,[1,7,41,287],2>, <291,[1,3,97,291],2>,
  <299,[1,13,23,299],2>, <303,[1,3,101,303],2>, <321,[1,3,107,321],2>, <339,[1,3,113,339],2>,
  <341,[1,11,31,341],2>, <371,[1,7,53,371],2>, <415,[1,5,83,415],2>, <417,[1,3,139,417],2>,
  <453,[1,3,151,453],2>, <471,[1,3,157,471],2>,
  <1155,[1,3,5,15,77,231,385,1155],2>, <1155,[1,7,11,15,77,105,165,1155],2>,
  <1155,[1,3,7,21,55,165,385,1155],2>,
  <1365,[1,3,5,7,13,15,21,35,39,65,91,105,195,273,455,1365],2>,
  <57,[1],3>, <145,[1,145],3>, <115,[1,5],4>, <143,[1,11],4>, <237,[1,237],4>, <93,[1],5>,
  <143,[1,13],6>, <161,[1,7],6>, <215,[1,43],6>, <237,[1,3],7>, <215,[1,5],8> ];

nin := 0;
for t in targets do
    D := t[1]; Ws := t[2]; g := t[3];
    if #Ws gt 2 then
        Say(Sprintf("D = %o W = %o g = %o: OUT -- W is bigger than <w_m>, so every Atkin-Lehner"
            * " involution is trivial on the quotient and none can be the hyperelliptic one", D, Ws, g));
        continue;
    end if;
    m := Maximum(Ws);
    wd := MolWeierstrassData(D, m, g);
    if #wd ne 1 then
        Say(Sprintf("D = %o W = %o g = %o: OUT -- no divisor n of D has 2g+2 = %o fixed points"
            * " on the quotient (%o candidates)", D, Ws, g, 2*g+2, #wd));
        continue;
    end if;
    n := wd[1][1]; np := wd[1][2]; rows := wd[1][3];
    if not (D in {m, n, np}) then
        Say(Sprintf("D = %o W = %o g = %o: OUT -- hyperelliptic involution w_%o, and D is not in"
            * " %o, so the points at infinity are not a conjugate pair", D, Ws, g, n, {m,n,np}));
        continue;
    end if;
    degs := Sort(&cat[[r[2] : i in [1..r[3]]] : r in rows]);
    reach := Maximum(degs) le 4;
    nin +:= reach select 1 else 0;
    Say(Sprintf("D = %o W = %o g = %o: IN (w_%o, w_%o), factor degrees %o -- %o",
        D, Ws, g, n, np, degs,
        reach select "within IndexFormEquation's reach"
              else Sprintf("BEYOND: a factor of degree %o", Maximum(degs))));
end for;
Say(Sprintf("IN SCOPE AND REACHABLE: %o of %o", nin, #targets));
Say("SCOPE DONE");
quit;
