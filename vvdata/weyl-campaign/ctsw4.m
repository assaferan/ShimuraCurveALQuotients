// GAP 5, step 5: the closed form for D odd, N odd -- tested PREDICTIVELY.
//
// Read off ctsw3.log (D odd, N odd; fitted on N = 3 and N = 7 ONLY):
//
//   -a_E(m) = H(4m) * prod_{p|D} (1 - (-m|p))/(p-1) * 12 (N - (-m|N))/(N^2 - 1)
//
// for m squarefree.  The Eichler factor came from c_3 = 3/2 and c_7 = 1/4, i.e.
// two points through the one-parameter guess 12/(N^2-1) -- which is exactly the
// kind of fit this campaign has been burned by.  So test it at N = 5, 11, 13,
// none of which was used, and at fresh D as well.
//
// Note the ramified factor (1 - (-m|p)) vanishes precisely when p SPLITS in
// Q(sqrt(-m)), so the support rule is inside the formula.
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/ctlat.m";
DEPTH := 120;

Hur := function(n)
    if n eq 0 then return Rationals()!(-1)/12; end if;
    if n mod 4 in {1, 2} then return Rationals()!0; end if;
    s := Rationals()!0; f := 1;
    while f^2 le n do
        if n mod f^2 eq 0 then
            d := -(n div f^2);
            if d mod 4 in {0, 1} then
                w := d eq -3 select 6 else (d eq -4 select 4 else 2);
                s +:= Rationals()!(2*ClassNumber(QuadraticForms(d)))/w;
            end if;
        end if;
        f +:= 1;
    end while;
    return s;
end function;
T3 := ThetaSeries(StandardLattice(3), 4*DEPTH);
error if exists{ n : n in [1..90] |
    (n mod 4 in {1,2} and Coefficient(T3,n) ne 12*Hur(4*n)) or
    (n mod 8 eq 3 and Coefficient(T3,n) ne 24*Hur(n)) }, "Hurwitz routine wrong";
printf "Hurwitz gate: passed\n\n";

massmul := function(Dp, Rl)
    ps := PrimeDivisors(Dp);
    m := Rationals()!(&*[ Integers() | p-1 : p in ps ]) / (48*2^(#ps-1));
    for p in PrimeDivisors(Rl) do m *:= Rationals()!(p+1)/2; end for;
    return m;
end function;
thetabar := function(Dp, Rl)
    A := ChangeRing(grossGram(Dp, Rl), Rationals());
    L := LatticeWithGram(A);
    reps := Representatives(Genus(L));
    raw := [ Rationals() | 0 : m in [0..2*DEPTH] ];
    den := Rationals()!0;
    for Lr in reps do
        w := Rationals()!1/#AutomorphismGroup(Lr);
        T := ThetaSeries(Lr, 2*DEPTH);
        for m in [0..2*DEPTH] do raw[m+1] +:= w*Coefficient(T, m); end for;
        den +:= w;
    end for;
    raw := [ x/den : x in raw ];
    return [ raw[2*m+1] : m in [0..DEPTH-1] ];
end function;

//        D,  N, D', R, Rs        N used in the fit?
cases := [ <15,7,5,7,21,  true>,  <35,3,7,3,15,  true>,     // the two fitted
           <21,5,7,5,15,  false>, <33,5,11,5,15, false>,    // N = 5  FRESH
           <15,11,5,11,33,false>,                           // N = 11 FRESH
           <15,13,5,13,39,false> ];                         // N = 13 FRESH

for cs in cases do
    D := cs[1]; N := cs[2]; Dp := cs[3]; Rl := cs[4]; Rls := cs[5];
    fitted := cs[6]; DN := D*N;
    m1 := massmul(Dp, Rl); msx := massmul(Dp, Rls);
    t1 := thetabar(Dp, Rl); ts := thetabar(Dp, Rls); t3 := thetabar(DN, 1);
    w1 := -m1/(2*(msx-m1)); ws := msx/(2*(msx-m1));
    G := [ Rationals() | ws*ts[i] + w1*t1[i] - (1/2)*t3[i] : i in [1..DEPTH] ];

    nok := 0; nbad := 0; bl := [ ];
    for m in [1..DEPTH-1] do
        if not IsSquarefree(m) then continue; end if;
        pred := -Hur(4*m)
                * (&*[ Rationals() | (1 - KroneckerSymbol(-m,p))/(p-1)
                       : p in PrimeDivisors(D) ])
                * (Rationals()!12*(N - KroneckerSymbol(-m,N))/(N^2-1));
        if G[m+1] eq pred then nok +:= 1;
        else nbad +:= 1; if #bl lt 4 then Append(~bl, <m, G[m+1], pred>); end if;
        end if;
    end for;
    printf "%o_%-3o  N = %-3o %-9o  %o of %o squarefree m %o\n",
        D, N, N, fitted select "(fitted)" else "(FRESH)",
        nok, nok+nbad,
        nbad eq 0 select "agree   PREDICTED"
        else Sprintf("agree, %o FAIL  e.g. %o", nbad, bl);
end for;
printf "DONE\n";
quit;
