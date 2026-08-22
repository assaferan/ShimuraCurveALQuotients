// Steps (ii) and the constancy character of the theta_g proof:
//  (A) the word unit u(w) = rho~(word)/rho_D(matrix) equals the Kubota cocycle
//      product along the word, computed EXACTLY via Stromberg Thm 4.1:
//         sigma(A, B) = mu(A, B) s(A) s(B) s(AB)^{-1},
//         mu(A, B) = (sigma_A sigma_AB, sigma_B sigma_AB)_infty,
//         sigma_M = c if c ne 0 else d;  s(M) = 1 if c ne 0 else sign(d).
//      rho~(word) = product of canonical generator lifts, so
//         (word matrix lift) = [prod of cocycles] * (canonical lift), i.e.
//         u(w) = prod_{prefix} sigma(prefix, next letter)   (as +-1).
//  (B) e_0 is a rho(Gamma_0(levD))-eigenvector: for gamma in Gamma_0(levD) with
//      c = levD * k, rho~(gamma-word) e_0 = chi(gamma) e_0 -- verified numerically
//      for a batch of gamma, and chi is recorded (the constancy character).
//
//   magma -b DD:=15 NN:=2 rho6.m
AttachSpec("ShimuraQuotients.spec");

D := 15; N := 2;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;

PREC := 60;
CC := ComplexField(PREC); ii := CC.1; pi := Pi(CC);

Ld := ShimuraCurveLattice(D, N);
M := IsOdd(D*N) select 4*D*N else 2*D*N;
fftdata := VVWeilFFT(Ld, CC : Dual := true);
elts := fftdata[7]; i0 := fftdata[8];
n := #elts;
printf "BASE %o %o  M = %o  |D| = %o\n", D, N, M, n;

// ---- (A) Kubota cocycle product per word ----
sigM := func< A | A[2][1] ne 0 select A[2][1] else A[2][2] >;
sfn := func< A | A[2][1] ne 0 select 1 else Sign(A[2][2]) >;
hilb := func< x, y | (x lt 0 and y lt 0) select -1 else 1 >;
sigma_cocycle := function(A, B)
    AB := A*B;
    mu := hilb(sigM(A)*sigM(AB), sigM(B)*sigM(AB));
    return mu * sfn(A) * sfn(B) * sfn(AB);
end function;

Smat := Matrix(Integers(), 2, 2, [0, -1, 1, 0]);
Tmat := Matrix(Integers(), 2, 2, [1, 1, 0, 1]);

// the letter matrices of a word, in application order (VVRhoInvE0FFT applies
// the INVERSE lifts letterwise for word = [l1, ..., lk] with matrix l1...lk)
// rho~(word)^{-1} = rho~(lk)^{-1} ... rho~(l1)^{-1}; the cocycle product relating
// the product of canonical lifts to the canonical lift of the product telescopes:
//   prod_i (canonical lift of l_i) = [prod of sigma(prefix_i, l_{i+1})] (canonical of prod)
kubota_word := function(word)
    mats := [ ];
    for t in word do
        if t[1] eq "S" then Append(~mats, Smat);
        else Append(~mats, Tmat^(t[2])); end if;
    end for;
    u := 1;
    P := mats[1];
    for i := 2 to #mats do
        u *:= sigma_cocycle(P, mats[i]);
        P := P*mats[i];
    end for;
    return u, P;
end function;

reps := VVCosetReps(M);
words := [ VVSTWord(g) : g in reps ];
// measured u from rho5-style comparison would need stromberg; here instead verify
// the INTERNAL consistency: rho~(word) computed letterwise (VVRhoInvE0FFT) vs
// rho~(word) predicted as kubota * rho of canonical lift -- we use the rho5-verified
// stromberg values only through the CLAIM u = kubota; so verify:
//   VVRhoInvE0FFT(word) = kubota(inverse word letters) * stromberg(matrix^{-1})
// To stay self-contained here, just print kubota units per coset class for the
// canonical reps and ALL cosets; cross-check against rho5's measured u happens in
// the analysis (the rho5 outputs record u per coset).
printf "KUBOTA units (word for w; unit for the INVERSE-letter word):\n";
for wi->w in words do
    gmat := VVWordMatrix(w);
    // VVRhoInvE0FFT applies inverse lifts of letters l1..lk in order l1..lk --
    // this equals rho~ of the word [lk^-1 ... l1^-1] read as a product; compute
    // kubota for that letter sequence:
    inv := [ ];
    for i := #w to 1 by -1 do
        t := w[i];
        if t[1] eq "S" then Append(~inv, <"S", 0>);   // S^{-1} = Z^2 S-ish: handle below
        else Append(~inv, <"T", -t[2]>); end if;
    end for;
    // S^{-1} = -S as a matrix; canonical lift bookkeeping for the inverse letters is
    // exactly what VVRhoInvE0FFT implements via conjugation; for the unit we only
    // need the matrix-level cocycles of the sequence of matrices applied:
    mats := [ ];
    for t in inv do
        if t[1] eq "S" then Append(~mats, -Smat);   // S^{-1}
        else Append(~mats, Tmat^(t[2])); end if;
    end for;
    u := 1; P := mats[1];
    for i := 2 to #mats do
        u *:= sigma_cocycle(P, mats[i]);
        P := P*mats[i];
    end for;
    cw := gmat[2][1] mod M; dw := gmat[2][2] mod M;
    printf "KUB c=%o d=%o gcd=%o unit=%o\n", cw, dw, GCD(cw, M), u;
end for;

// ---- (B) the constancy character: gamma in Gamma_0(levD) acts on e_0 by chi ----
Qr := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
frac := func< r | r - Floor(r) >;
lev := 1;
Qbar := [ frac(-((ChangeRing(elts[i]@@Ld`to_disc, Rationals()))*Qr,
                 ChangeRing(elts[i]@@Ld`to_disc, Rationals()))/(2*dn^2)) : i in [1..n] ];
for L in Divisors(2*M^2) do
    if forall{ j : j in [1..n] | frac(L*Qbar[j]) eq 0 } then lev := L; break; end if;
end for;
printf "level = %o\n", lev;
printf "CHI values on Gamma_0(%o) (gamma = [a b; lev*k d]):\n", lev;
cnt := 0;
for k in [0, 1, -1, 2] do
    for d0 in [1, -1, 5, 7, 11, 13] do
        c := lev*k;
        // solve a d - b c = 1 with small entries
        founda := false;
        for a in [-30..30] do
            if c eq 0 then
                if a*d0 eq 1 then bb := 0; founda := true; end if;
            else
                if (a*d0 - 1) mod c eq 0 then bb := (a*d0 - 1) div c; founda := true; end if;
            end if;
            if founda then
                gam := Matrix(Integers(), 2, 2, [a, bb, c, d0]);
                wgam := VVSTWord(gam);
                rv := VVRhoInvE0FFT(fftdata, wgam);   // rho~(word for gam)^{-1}-style e_0
                // eigen check: all mass on e_0?
                mx := Maximum([ Abs(rv[i]) : i in [1..n] | i ne i0 ]);
                chi := rv[i0];
                printf "  gamma=(%o,%o;%o,%o): offdiag %o  chi=(%o,%o) arg8=%o\n",
                    a, bb, c, d0, RealField(4)!mx,
                    RealField(8)!Re(chi), RealField(8)!Im(chi),
                    RealField(8)!(8*Arg(chi)/(2*Pi(RealField(20))));
                cnt +:= 1;
                break;
            end if;
        end for;
        if cnt ge 16 then break; end if;
    end for;
    if cnt ge 16 then break; end if;
end for;
quit;
