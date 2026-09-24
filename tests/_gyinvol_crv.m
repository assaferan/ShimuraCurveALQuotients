// tests/_gyinvol_crv.m -- NOT a test (leading underscore: excluded from the CI matrix). The CRV
// counterpart of tests/_gyinvol.m; run by hand to REGENERATE the 4x4 ws_data matrices in
// tests/X0_{26_3,57_1,21_2,14_3}.m. Committed as their PROVENANCE.
//
// Transport Guo-Yang's published involutions onto our CRV (paired) models at W={1}.
// Same logic as tests/_gyinvol.m, but psi comes from construct_crv_isomorphism -- IsIsomorphic
// HANGS on paired presentations (>1 h at 26_3).
AttachSpec("ShimuraQuotients.spec");
import "tests/_crviso.m" : construct_crv_isomorphism;

// ⚠ IdentityMap on a WEIGHTED ambient is a TorMap, so `phi eq IdentityMap(C)` raises "Bad argument
// types". Check the two properties algebraically instead.
// (a) M preserves C: every defining polynomial, substituted, lies in C's ideal.
function preserves(C, M)
    R := CoordinateRing(AmbientSpace(C));
    img := Eltseq(Vector([R.i : i in [1..Rank(R)]])*ChangeRing(M, R));
    I := Ideal(C);
    return &and[ NormalForm(Evaluate(f, img), I) eq 0 : f in DefiningPolynomials(C) ];
end function;
// (b) M is an involution PROJECTIVELY: M^2 must act as (x_i) -> (lambda^{w_i} x_i) for one lambda,
// not as the literal identity -- on P(1,3,1,1), diag(-1,-1,-1,1) squares to the identity but
// diag(-1,-1,-1,1) itself is already the identity map projectively.
function is_involution(C, M, wts)
    N := M^2;
    for i, j in [1..Nrows(N)] do
        if (i ne j) and (N[i][j] ne 0) then return false; end if;
    end for;
    lam := 0;
    for i in [1..Nrows(N)] do if wts[i] eq 1 then lam := N[i][i]; break; end if; end for;
    if lam eq 0 then return false; end if;
    return &and[ N[i][i] eq lam^wts[i] : i in [1..Nrows(N)] ];
end function;

// base, GY weights, GY equations as f(a,b,c,d), [<label, 4x4 matrix on (a,b,c,d)>]
data := [*
  // a=x, b=y (wt 2), c=z, d=hom.   z^2=-9x^2-2, y^2=-7x^4+22x^2+1
  <"14_3", [1,2,1,1],
   func<a,b,c,d | [c^2 + 9*a^2 + 2*d^2, b^2 + 7*a^4 - 22*a^2*d^2 - d^4]>,
   [* <2,  DiagonalMatrix(Rationals(),[-1, 1, 1,1])>,   // (-x, y, z)
      <3,  DiagonalMatrix(Rationals(),[ 1,-1,-1,1])>,   // ( x,-y,-z)
      <14, DiagonalMatrix(Rationals(),[ 1,-1, 1,1])> *]>,   // ( x,-y, z)
  // a=x, b=y (wt 3), c=z, d=hom.   z^2=-8x^2-3, y^2=x^6-2x^4+9x^2+8
  <"26_3", [1,3,1,1],
   func<a,b,c,d | [c^2 + 8*a^2 + 3*d^2, b^2 - (a^6 - 2*a^4*d^2 + 9*a^2*d^4 + 8*d^6)]>,
   [* <2,  DiagonalMatrix(Rationals(),[-1,-1,-1,1])>,   // (-x,-y,-z)
      <3,  DiagonalMatrix(Rationals(),[ 1,-1,-1,1])>,   // ( x,-y,-z)
      <26, DiagonalMatrix(Rationals(),[ 1,-1, 1,1])> *]>,   // ( x,-y, z)
  // a=x, b=y (wt 3), c=z, d=hom.   z^2=-x^2-3, y^2=-(3x-1)(3x+1)(x^2+7)(x^2+3)
  <"21_2", [1,3,1,1],
   func<a,b,c,d | [c^2 + a^2 + 3*d^2,
                   b^2 + (3*a-d)*(3*a+d)*(a^2+7*d^2)*(a^2+3*d^2)]>,
   [* <2, DiagonalMatrix(Rationals(),[-1,-1,-1,1])>,    // (-x,-y,-z)
      <3, DiagonalMatrix(Rationals(),[ 1, 1,-1,1])>,    // ( x, y,-z)
      <7, DiagonalMatrix(Rationals(),[ 1,-1, 1,1])> *]>,    // ( x,-y, z)
  // ⚠ GY's variables here are (s,x,y): a=s, b=y (wt 2), c=x, d=hom.
  //    y^2=(3s+1)(3s^3+11s^2+17s+1), x^2=-4s^2+2s-1
  <"57_1", [1,2,1,1],
   func<a,b,c,d | [b^2 - (3*a+d)*(3*a^3+11*a^2*d+17*a*d^2+d^3), c^2 + 4*a^2 - 2*a*d + d^2]>,
   [* <19, DiagonalMatrix(Rationals(),[1,-1, 1,1])>,    // (s, x,-y)
      <57, DiagonalMatrix(Rationals(),[1, 1,-1,1])> *]> // (s,-x, y)
*];

for d0 in data do
    b, wgy, egy, invs := Explode(d0);
    printf "\n================ %o ================\n", b;
    models := eval (Read("data/models/models_" cat b cat ".m") cat "\nreturn models;");
    e := models[[Integers()|1]][1];
    gour := e[1]; strs := e[3];
    R<yy,xx,ss,zz> := PolynomialRing(Rationals(), 4);
    dy := Degree(eval ("return " cat strs[1] cat ";"))
          where y is yy where x is xx where s is ss where z is zz;
    Po<xo,yo,so,zo> := WeightedProjectiveSpace(Rationals(), [1, dy div 2, 1, 1]);
    Cours := Curve(Po, [eval ("return " cat st cat ";") : st in strs]
                   where y is yo where x is xo where s is so where z is zo);
    Qg<a,b_,c,dd> := WeightedProjectiveSpace(Rationals(), wgy);
    Cgy := Curve(Qg, egy(a,b_,c,dd));
    printf "  our genus %o (y wt %o), GY genus %o (wts %o)\n",
           Genus(Cours), dy div 2, Genus(Cgy), wgy;
    okp, psi := construct_crv_isomorphism(Cours, Cgy);
    if not okp then
        // ⚠ The constructor declines when the two pairs present the curve over DIFFERENT
        // intermediate quotients (21_2: GY's y has weight 3, a genus-2 y-quotient, ours weight 2
        // and genus 1), so it cannot take a Mobius map from a common base. Fall back to the
        // general IsIsomorphic -- correct but slow on paired presentations, which is why the
        // constructor exists at all.
        printf "  constructor declined; falling back to IsIsomorphic (slow)...\n";
        t0 := Cputime();
        okp, psi := IsIsomorphic(Cours, Cgy);
        printf "  IsIsomorphic: %o after %o s\n", okp, Cputime(t0);
        if not okp then printf "  NOT ISOMORPHIC\n"; continue; end if;
    end if;
    Rg := CoordinateRing(Qg); cg := [Rg.i : i in [1..4]];
    Ro := CoordinateRing(Po);  co := [Ro.i : i in [1..4]];
    for iv in invs do
        Q, M := Explode(iv);
        wgy_map := map< Cgy -> Cgy | Eltseq(Vector(cg)*ChangeRing(M, Rg)) >;
        printf "  w_%o: preserves GY's curve? %o; involution there? %o\n",
               Q, preserves(Cgy, M), is_involution(Cgy, M, wgy);
        // ⚠ Inverse(psi) raises "Map has no inverse" here even though psi IS an isomorphism --
        // Magma cannot invert this particular representation. Try the routes that do work.
        if not assigned psinv then
            for route in ["IsInvertible", "Expand", "IsIsomorphism"] do
                try
                    case route:
                        when "IsInvertible":  okv, cand := IsInvertible(psi);
                        when "Expand":        cand := Inverse(Expand(psi)); okv := true;
                        when "IsIsomorphism": okv, cand := IsIsomorphism(psi);
                    end case;
                    if okv then psinv := cand; printf "  inverse via %o\n", route; break; end if;
                catch e printf "  %o failed: %o\n", route, e`Object; end try;
            end for;
        end if;
        if not assigned psinv then printf "  NO INVERSE BY ANY ROUTE\n"; break; end if;
        ourw := psi*wgy_map*psinv;
        // ⚠ SOLVE for the matrix rather than reading it off. The composite is a single
        // unreduced degree-39 representation, so Exponents() sees a non-linear form even when the
        // MAP is linear. On P(1,wy,1,1) only y has weight wy, so a weight-respecting matrix must
        // send y -> c*y and act on (x,s,z) by a 3x3 block. Write the composite's images as
        // q1,q2,q3,q4; the matrix map agrees with it iff there is lambda with q1 = lambda*L1,
        // q3 = lambda*L3, q4 = lambda*L4 for LINEAR L's, i.e. iff
        //     q1*L3 - q3*L1  and  q1*L4 - q4*L1  vanish on the curve.
        // Those are LINEAR conditions on the 9 unknown coefficients, so the L's are a kernel.
        de0 := DefiningEquations(ourw);
        I := Ideal(Cours);
        wt1 := [co[1], co[3], co[4]];                       // x, s, z -- the weight-1 coordinates
        N1 := [NormalForm(de0[1]*m, I) : m in wt1];
        N3 := [NormalForm(de0[3]*m, I) : m in wt1];
        N4 := [NormalForm(de0[4]*m, I) : m in wt1];
        // unknowns (a1,a2,a3 | b1,b2,b3 | c1,c2,c3) for (L1 | L3 | L4)
        rows := []; monset := {};
        for pair in [<N3, N1, 2>, <N4, N1, 3>] do
            Na, Nb, blk := Explode(pair);
            for i in [1..3] do monset join:= Set(Monomials(Na[i])) join Set(Monomials(Nb[i])); end for;
        end for;
        monlist := Setseq(monset);
        eqs := [];
        for pair in [<N3, N1, 2>, <N4, N1, 3>] do
            Na, Nb, blk := Explode(pair);
            for mm in monlist do
                row := [Rationals()| 0 : k in [1..9]];
                for i in [1..3] do
                    row[i]       -:= MonomialCoefficient(Na[i], mm);   // -a_i * q3(or q4)*mon_i
                    row[3*(blk-1)+i] +:= MonomialCoefficient(Nb[i], mm);
                end for;
                Append(~eqs, row);
            end for;
        end for;
        K := Kernel(Transpose(Matrix(Rationals(), eqs)));
        printf "    linear solve: %o monomial equation(s), kernel dimension %o\n", #eqs, Dimension(K);
        if Dimension(K) eq 0 then printf "    NO LINEAR MAP EXISTS -- ws_data cannot hold it\n"; continue; end if;
        v := Eltseq(K.1);
        Mo := ZeroMatrix(Rationals(), 4, 4);
        idx := [1,3,4];
        for j in [1..3] do for i in [1..3] do Mo[idx[i]][idx[j]] := v[3*(j-1)+i]; end for; end for;
        // y -> c*y: only the sign can differ once the block is fixed, so test both against the MAP
        got := false;
        for cc in [1,-1] do
            Mo[2][2] := cc;
            try
                chk := map< Cours -> Cours | Eltseq(Vector(co)*ChangeRing(Mo, Ro)) >;
                if chk eq ourw then got := true; break; end if;
            catch e ; end try;
        end for;
        if not got then printf "    solved block does NOT reproduce the map\n"; continue; end if;
        printf "    ws_data[{1}][%o] := Matrix(4,4,%o);\n", Q,
               [Rationals()| Mo[i][j] : j in [1..4], i in [1..4]];
        printf "        // ==the transported MAP: %o; preserves ours: %o; involution: %o\n",
               chk eq ourw, preserves(Cours, Mo),
               is_involution(Cours, Mo, [1, dy div 2, 1, 1]);
    end for;
end for;
exit;
