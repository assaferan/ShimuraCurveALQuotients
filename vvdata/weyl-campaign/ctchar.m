// LEMMA, half 2, done properly: MEASURE the character of each Gross theta.
//
// Pool membership (ctpoolgen.log) confirms the character wherever it succeeds,
// but the eta pool is only a SUBSPACE of M_{3/2}, so non-membership is
// inconclusive -- 12 of 26 instances came back inconclusive that way.  This
// computes the character directly instead.
//
// For gamma in Gamma_0(l), theta_L | gamma = chi(gamma) theta_L, and a_0(theta_L)
// = 1, so a_0(theta_L | gamma) IS chi(gamma).  ctTheta evaluates exactly that
// closed form, so feeding it words for gamma in Gamma_0(l) reads the character
// off directly.
//
// PREDICTION worth stating in advance: det A = 2(Delta R)^2 (measured, 100/100
// in ctlevel.log), so 2 det A = (2 Delta R)^2 is a perfect square.  If the
// classical character of a ternary theta is a Kronecker symbol attached to
// 2 det A, it should therefore collapse to something depending only on d mod 8
// -- the SAME character for every (Delta, R).  Candidates tested: 1, (-1/d),
// (2/d), (-2/d).
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/fastcosets.m";
load "vvdata/weyl-campaign/etared.m";
load "vvdata/weyl-campaign/ctlat.m";
PREC := 60; CC := ComplexField(PREC); QSIGN := -1;

oddpart := function(n) m := n; while IsEven(m) do m div:= 2; end while; return m; end function;

// (Delta, R) drawn from bases that PASSED pool membership and from those that
// did not, so the two groups can be compared directly.
cases := [ <5,2>, <5,6>, <30,1>,        // 15_2   (passed)
           <7,2>, <7,6>, <42,1>,        // 21_2   (inconclusive)
           <11,3>, <11,6>, <66,1>,      // 22_3   (inconclusive)
           <11,2>, <66,1>,              // 33_2   (passed)
           <11,10>, <110,1>,            // 55_2   (inconclusive)
           <5,7>, <5,14>, <70,1>,       // 10_7   (passed)
           <11,1>, <13,1>, <7,1> ];     // N = 1 lattices

printf "%-12o %-8o %-6o  %o\n", "lattice", "level", "d", "a_0(theta|gamma)  vs  1 / (-1|d) / (2|d) / (-2|d)";
for c in cases do
    Dl := c[1]; Rl := c[2];
    if #PrimeDivisors(Dl) mod 2 eq 0 then continue; end if;
    A := ChangeRing(grossGram(Dl, Rl), Rationals());
    lev := 4*oddpart(Dl*Rl);
    verdict := [ true, true, true, true ];   // 1, (-1|d), (2|d), (-2|d)
    for d in [1..40] do
        if GCD(d, 2*lev) ne 1 then continue; end if;
        // gamma = [a b; lev d] in Gamma_0(lev) with a*d - b*lev = 1
        a := InverseMod(d mod lev, lev);
        if (a*d - 1) mod lev ne 0 then continue; end if;
        b := (a*d - 1) div lev;
        gam := Matrix(Integers(), 2, 2, [a, b, lev, d]);
        if Determinant(gam) ne 1 then continue; end if;
        val := ctTheta(A, [ VVSTWord(gam) ], CC, QSIGN)[1];
        cands := [ CC | 1, KroneckerSymbol(-1, d), KroneckerSymbol(2, d),
                        KroneckerSymbol(-2, d) ];
        for j in [1..4] do
            if Abs(val - cands[j]) gt 10^(-30) then verdict[j] := false; end if;
        end for;
        if d le 13 then
            printf "%-12o %-8o %-6o  %o\n", Sprintf("Gr(%o,%o)", Dl, Rl), lev, d,
                ComplexField(8)!val;
        end if;
    end for;
    names := [ "1", "(-1|d)", "(2|d)", "(-2|d)" ];
    printf "   --> Gr(%o,%o) level %o: character is %o\n", Dl, Rl, lev,
        #[ j : j in [1..4] | verdict[j] ] eq 0 select "NONE OF THE FOUR"
        else [ names[j] : j in [1..4] | verdict[j] ];
end for;
printf "DONE\n";
quit;
