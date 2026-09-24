// The two-coset local Whittaker density, efficient version: all-integer counting, character handled
// by accumulating counts per residue class of the phase (one small cyclotomic combination at the
// end).  This file only DEFINES the function and self-validates it on X0^15(2); the assembly probe
// loads it with "load".
//
// W_{m,p}(mu; eta*) = lim_k p^{-2k} Sum_{y in (eta*+L)/p^kL, Q(y) = m mod p^k} psi(B(mu, y))
// with Q(y) = y G y^t / 2, B(u,v) = u G v^t, G = the Gram the pipeline uses (negQ), psi = e(.).

function PPart(x, p, slack)
    // the class of x mod Z_p, as a rational with p-power denominator
    b := Denominator(x); a := Numerator(x);
    vp := Valuation(b, p);
    if vp eq 0 then return Rationals()!0; end if;
    bp := b div p^vp;
    return Rationals()!((a * InverseMod(bp, p^(vp+slack))) mod p^(vp+slack)) / p^vp;
end function;

function TwoCosetW(m, p, k, seed, twist, G)
    // p-normalize m: its prime-to-p denominator is a p-adic unit; invert it so that congruences
    // mod p^k make sense.  (Cusp-0 indices r = j/M carry such denominators.)
    bm := Denominator(m); am := Numerator(m);
    vm := Valuation(bm, p); bmp := bm div p^vm;
    if bmp ne 1 then
        m := Rationals()!(am * InverseMod(bmp, p^(k + vm + 10))) / p^vm;
    end if;
    // p-reduce both cosets; scale the seed's denominator away so everything is integral
    sv := Vector(Rationals(), [PPart(x, p, k+4) : x in Eltseq(seed)]);
    tv := Vector(Rationals(), [PPart(x, p, k+4) : x in Eltseq(twist)]);
    e := Maximum([Valuation(Denominator(x), p) : x in Eltseq(sv)] cat [0]);
    Gq := ChangeRing(G, Rationals());
    // phases B(tv, y) for y in seed + Z^3 lie in (1/p^f)Z: bound f
    f := Maximum([Valuation(Denominator((tv*Gq, Vector(Rationals(), [i eq j select 1 else 0 : i in [1..3]]))), p)
                  : j in [1..3]] cat [Valuation(Denominator((tv*Gq, sv)), p)] cat [0]);
    ord := p^f;
    K := f eq 0 select Rationals() else CyclotomicField(ord);
    counts := [0 : i in [1..ord]];              // counts[r+1] for phase class r/ord
    pk := p^k; mod2 := p^(k + 2*e);
    // integer seed: S = p^e * sv
    S := [Integers()!(p^e * x) : x in Eltseq(sv)];
    // require m's p-denominator handled: M2 = p^{2e} * m must be p-integral for solutions to exist
    m2 := p^(2*e) * m;
    if Denominator(2*m2) ne 1 then return K!0, 0; end if;   // deeper p-denominator than the coset: no solutions
    // Gram doubled to avoid /2: Q(y) = (y G y)/2; with Y = p^e y: Y G Y = 2 p^{2e} m mod 2 p^{k+2e}
    G2 := ChangeRing(G, Integers());
    tGrow := [ (tv*Gq)[j] : j in [1..3] ];      // row vector t G
    tshift := (tv*Gq, sv);                       // constant part of the phase from the seed
    total := 0;
    for y1 in [0..pk-1], y2 in [0..pk-1], y3 in [0..pk-1] do
        Y := [ S[1] + p^e*y1, S[2] + p^e*y2, S[3] + p^e*y3 ];
        q2 := &+[ Y[i] * &+[ G2[i][j]*Y[j] : j in [1..3] ] : i in [1..3] ];   // = 2 p^{2e} Q(y)
        if (q2 - Integers()!(2*m2)) mod (2*mod2) eq 0 or (q2 - 2*m2) eq 0 then
            if f eq 0 then
                total +:= 1;
            else
                ph := tshift + tGrow[1]*y1 + tGrow[2]*y2 + tGrow[3]*y3;
                r := Integers()!( (Numerator(ph) * (ord div Denominator(ph))) mod ord );
                counts[r+1] +:= 1;
            end if;
        end if;
    end for;
    if f eq 0 then return K!total / p^(2*k), total; end if;
    z := K.1;
    return (&+[ counts[r+1] * z^r : r in [0..ord-1] ]) / p^(2*k), &+counts;
end function;

// ---- self-validation when run standalone -------------------------------------------------------
if assigned selfcheck then
    D := 15; N := 2;
    Ld := ShimuraCurveLattice(D, N);
    negQ := -ChangeRing(Ld`Q, Integers());
    Lfull := RSpaceWithBasis(IdentityMatrix(Integers(), 3));
    denom := Ld`denom; to_disc := Ld`to_disc;
    G := Ld`disc_grp;
    zero := Vector(Rationals(), [0,0,0]);
    all := Sort([ Eltseq(ChangeRing(g@@to_disc, Rationals())/denom) : g in G ]);
    picks := [ zero ];
    for p in [2, 3, 5] do
        got := 0;
        for v in all do
            if got ge 2 then break; end if;
            if Lcm([Denominator(x) : x in v] cat [1]) mod p eq 0 then
                w := Vector(Rationals(), v);
                if not w in picks then Append(~picks, w); got +:= 1; end if;
            end if;
        end for;
    end for;
    nbad := 0; nok := 0;
    for w in picks do
        for p in [2, 3, 5] do
            kk := p eq 2 select 5 else (p eq 3 select 4 else 3);
            for m in [1..6] do
                code := LocalWhittakerAtOne(Rationals()!m, p, w, Lfull, negQ);
                dens := TwoCosetW(Rationals()!m, p, kk, w, zero, negQ);
                ok := IsCoercible(Rationals(), dens) and Rationals()!dens eq code;
                if ok then nok +:= 1; else nbad +:= 1;
                    printf "  DISAGREES w=%o p=%o m=%o: count=%o code=%o\n", w, p, m, dens, code;
                end if;
            end for;
        end for;
    end for;
    printf "seed-coset validation: %o ok, %o disagree\n", nok, nbad;
    // stability: k vs k+1 on a spot
    w := picks[2];
    for p in [2,3] do
        a := TwoCosetW(Rationals()!3, p, 3, w, zero, negQ);
        b := TwoCosetW(Rationals()!3, p, 4, w, zero, negQ);
        printf "stability p=%o m=3: k=3 -> %o, k=4 -> %o %o\n", p, a, b, a eq b select "stable" else "*** UNSTABLE";
    end for;
    quit;
end if;
