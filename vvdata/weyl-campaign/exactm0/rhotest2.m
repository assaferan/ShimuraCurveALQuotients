// Which factor of the exact per-class product fails for a representative with a != d?
// gamma = [3 1; 5 2] vs [1 0; 5 1] at X_0^15(2):
//   (A) rho row:   v_fft[eta] (ST-word lift) against e(1/8) c^{-3/2} |D|^{-1/2} e(d Q/c) G(c,a;eta), ALL eta,
//                  and against the a <-> d swapped version;
//   (B) slash constant per monomial: (c tau + d)^{-1/2} prod eta(d gamma tau)^{r_d} / sfun(tau) against
//                  e(-1/8) prod_d [eps(g_d) e(b_d/(24 e_d)) e_d^{-1/2}]^{r_d}, at two points tau.
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
D := 15; N := 2;
Ld := ShimuraCurveLattice(D, N);
CC := ComplexField(60); ii := CC.1; ee := func< z | Exp(2*Pi(CC)*ii*z) >;
fftdata := VVWeilFFT(Ld, CC : Dual := true); elts := fftdata[7]; i0 := fftdata[8];
Qr := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
vs := [ChangeRing(g@@Ld`to_disc, Rationals()) : g in elts];
nG := #elts;
Qe := [ (vs[i]*Qr, vs[i])/(2*dn^2) : i in [1..nG] ];

function gaussnum(c, a, i)
    G := CC!0;
    for nu in CartesianPower([0..c-1], 3) do
        nv := ChangeRing(Vector(Integers(), [nu[1], nu[2], nu[3]]), Rationals());
        G +:= ee((a*(nv*Qr, nv)/2 + (vs[i]*Qr, nv)/dn)/c);
    end for;
    return G;
end function;

for mm in [[1,0,5,1], [3,1,5,2], [2,1,3,2], [1,2,4,9], [7,2,10,3]] do
    g := Matrix(Integers(), 2, 2, mm);
    a := g[1][1]; c := g[2][1]; d := g[2][2];
    v := VVRhoInvE0FFT(fftdata, VVSTWord(g));
    pre := ee(1/8) * CC!c^(-3/2) / (D*N*Sqrt(CC!2));
    fd := [ pre * ee(d*Qe[i]/c) * gaussnum(c, a, i) : i in [1..nG] ];
    fa := [ pre * ee(a*Qe[i]/c) * gaussnum(c, d, i) : i in [1..nG] ];
    k := rep{i : i in [1..nG] | Abs(fd[i]) gt 10^-10};
    s := Round(Re(v[k]/fd[k]));
    errd := Maximum([Abs(v[i] - s*fd[i]) : i in [1..nG]]);
    k2 := rep{i : i in [1..nG] | Abs(fa[i]) gt 10^-10};
    s2 := Round(Re(v[k2]/fa[k2]));
    erra := Maximum([Abs(v[i] - s2*fa[i]) : i in [1..nG]]);
    printf "A %o: sign %o max|fft - formula(d, G(a))| = %o ; swapped (a, G(d)): sign %o err %o\n",
           mm, s, RealField(5)!errd, s2, RealField(5)!erra;
end for;

// (B) the slash constant
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fsa := BorcherdsForms(star, curves : Prec := 100);
f := fsa[-2]; R := Parent(f); ds := R`ds; M := R`M;
monos := Exponents(f);

function dedsum(h, k)
    s := Rationals()!0;
    for i := 1 to k - 1 do
        x := Rationals()!i/k; y := Rationals()!(h*i)/k; y := y - Floor(y);
        if y ne 0 then s +:= (x - 1/2)*(y - 1/2); end if;
    end for;
    return s;
end function;
function epsexp(g)
    a := g[1][1]; c := g[2][1]; d := g[2][2];
    x := (Rationals()!(a + d)/(12*c) + dedsum(-d, c))/2;
    return x - Floor(x);
end function;
function triang(g, d)
    c := g[2][1];
    g2 := Matrix(Integers(), 2, 2, [d*g[1][1], d*g[1][2], c, g[2][2]]);
    h := GCD(g2[1][1], c);
    p1 := g2[1][1] div h; p2 := c div h;
    _, u, v := XGCD(p1, p2);
    gd := Matrix(Integers(), 2, 2, [p1, -v, p2, u]);
    sd := gd^-1 * g2;
    return sd[1][1], sd[1][2], sd[2][2], gd;
end function;

for mm in [[1,0,5,1], [3,1,5,2], [2,1,3,2]] do
    g := Matrix(Integers(), 2, 2, mm);
    a := g[1][1]; b := g[1][2]; c := g[2][1]; d := g[2][2];
    tri := []; gds := [];
    for dd in ds do ad, bd, ed, gd := triang(g, dd); Append(~tri, <ad, bd, ed>); Append(~gds, gd); end for;
    nbad := 0; ntot := 0;
    for r in monos do
        kex := ee(-1/8) * &*[ CC | (ee(epsexp(gds[i])) * ee(Rationals()!tri[i][2]/(24*tri[i][3])) * CC!tri[i][3]^(-1/2))^r[i] : i in [1..#ds] ];
        for tau in [CC!0.31 + CC!1.31*ii, CC!(-0.57) + CC!1.73*ii] do
            gtau := (a*tau + b)/(c*tau + d);
            num := (c*tau + d)^(-1/2) * &*[ CC | DedekindEta(dd*gtau)^r[i] : i->dd in ds | r[i] ne 0 ];
            ord := &+[ Rationals() | r[i]*tri[i][1]/(24*tri[i][3]) : i in [1..#ds] ];
            sfun := ee(tau*ord) * &*[ CC | ( DedekindEta((tri[i][1]*tau + tri[i][2])/tri[i][3]) * ee(-(tri[i][1]*tau + tri[i][2])/(24*tri[i][3])) )^r[i] : i in [1..#ds] | r[i] ne 0 ];
            knum := num/sfun;
            ntot +:= 1;
            if Abs(knum - kex) gt 10^-30 * Maximum(1, Abs(kex)) then
                nbad +:= 1;
                if nbad le 3 then printf "  B %o r = %o: ratio num/exact = %o\n", mm, r, ComplexField(8)!(knum/kex); end if;
            end if;
        end for;
    end for;
    printf "B %o: slash constant exact vs numeric: %o bad of %o\n", mm, nbad, ntot;
end for;
quit;
