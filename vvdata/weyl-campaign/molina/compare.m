// Compare the equations the implementation produced with the ten Molina publishes.
// The comparison is by isomorphism over Q (a hyperelliptic model is only defined up to a change
// of coordinate) -- never by coefficients.
//   magma -b res:=<file of "D m  y^2 = ..." lines> logfile:=... compare.m

AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
P<x> := PolynomialRing(Rationals());

L := Open(logfile, "w");
procedure Say(s) fprintf L, "%o\n", s; Flush(L); end procedure;

molina := AssociativeArray();
molina[<39,1>]  := -(7*x^4+79*x^3+311*x^2+497*x+277)*(x^4+9*x^3+29*x^2+39*x+19);
molina[<55,1>]  := -(3*x^4-32*x^3+130*x^2-237*x+163)*(x^4-8*x^3+34*x^2-83*x+81);
molina[<35,5>]  := -x*(9*x+4)*(4*x+1)*(172*x^3+176*x^2+60*x+7);
molina[<51,17>] := -x*(7*x^3+52*x^2+116*x+68)*(x-1)*(x+3);
molina[<57,3>]  := -(x-9)*(x^3-19*x^2+119*x-249)*(7*x^2-104*x+388);
molina[<65,13>] := -(x^2-3*x+1)*(7*x^4-3*x^3-32*x^2+25*x-5);
molina[<65,5>]  := -(x^2+7*x+9)*(7*x^4+81*x^3+319*x^2+508*x+268);
molina[<69,23>] := -x*(x+4)*(4*x^4-16*x^3+11*x^2+10*x+3);
molina[<85,5>]  := -(3*x^2-41*x+133)*(x^4-23*x^3+183*x^2-556*x+412);
molina[<85,85>] := (x^2-3*x+1)*(x^4+x^3-15*x^2+20*x-8);

ours := eval Read(res);      // list of <D, m, [polynomials]>

nmatch := 0;
for e in ours do
    D := e[1]; m := e[2]; polys := e[3];
    key := <D, m>;
    C0 := HyperellipticCurve(molina[key]);
    hit := false;
    for f in polys do
        C := HyperellipticCurve(f);
        if IsIsomorphic(C, C0) then hit := true; break; end if;
    end for;
    // are all the equations we produced the same curve?
    same := &and[IsIsomorphic(HyperellipticCurve(f), HyperellipticCurve(polys[1])) : f in polys];
    Say(Sprintf("D = %o, w_%o: %o equation(s); isomorphic to Molina's: %o; all ours agree: %o",
        D, m, #polys, hit, same));
    if hit then nmatch +:= 1; end if;
end for;
Say(Sprintf("MATCHED %o of %o targets attempted (out of Molina's 10)", nmatch, #ours));
Say("COMPARE DONE");
quit;
