// Isolate the rho-formula failure of exactm0.m RHOCHECK at X_0^15(2):
//   (i)   v_fft      = rho*(gamma^-1) e_0 from VVRhoInvE0FFT (ST-word lift)
//   (ii)  v_grp      = row 0 of rho_B([1 0;1 1])^2, with rho_B([1 0;1 1])_{eta,beta} = const e(Q(eta+beta))
//                      (the c = 1 case of the theta-transformation formula, confirmed against S T^k S^-1 ...)
//   (iii) v_formula  = e(1/8) c^{-3/2} |D|^{-1/2} e(d Q(eta)/c) sum_{nu in L/cL} e((a Q(nu) + (eta,nu))/c)
// for gamma = [1 0; 2 1] = [1 0; 1 1]^2.  Ratios (i)/(ii) and (iii)/(ii) over a sample of eta.
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
D := 15; N := 2;
Ld := ShimuraCurveLattice(D, N);
CC := ComplexField(30); ii := CC.1; ee := func< z | Exp(2*Pi(CC)*ii*z) >;
fftdata := VVWeilFFT(Ld, CC : Dual := true); elts := fftdata[7]; i0 := fftdata[8];
Qr := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
vs := [ChangeRing(g@@Ld`to_disc, Rationals()) : g in elts];
nG := #elts;
Qe := [ (vs[i]*Qr, vs[i])/(2*dn^2) : i in [1..nG] ];
pair := func< i, j | (vs[i]*Qr, vs[j])/dn^2 >;
// index of eta + beta in elts
G := Ld`disc_grp; idxof := AssociativeArray(); for i->g in elts do idxof[g] := i; end for;
addidx := func< i, j | idxof[elts[i] + elts[j]] >;

g1 := Matrix(Integers(), 2, 2, [1, 0, 1, 1]);
g2 := g1^2;
v1 := VVRhoInvE0FFT(fftdata, VVSTWord(g1));
v2 := VVRhoInvE0FFT(fftdata, VVSTWord(g2));
// (ii): rho_B(g1)_{eta,beta} = K e(Q(eta + beta)); row 0 of g1: K e(Q(eta)); fix K from v1 (lift sign absorbed)
K := v1[i0];                                   // = K e(Q(0)) = K
printf "c = 1: v_fft/(K e(Q(eta))) over all eta: %o distinct values\n",
       #{ <Round(10^8*Re(z)), Round(10^8*Im(z))> where z := v1[i]/(K*ee(Qe[i])) : i in [1..nG] };
// row 0 of rho_B(g1)^2: sum_beta rho_{0,beta} rho_{beta,eta} = K^2 sum_beta e(Q(beta)) e(Q(beta + eta))
idx := [1 + (k*97) mod nG : k in [0..23]];
rnd := func< z | <Round(10^6*Re(z))/10^6, Round(10^6*Im(z))/10^6> >;
r12 := {}; r32 := {};
for i in idx do
    vg := K^2 * &+[ ee(Qe[b]) * ee(Qe[addidx(b, i)]) : b in [1..nG] ];
    a := 1; c := 2; d := 1;
    Gs := CC!0;
    for nu in CartesianPower([0..c-1], 3) do
        nv := ChangeRing(Vector(Integers(), [nu[1], nu[2], nu[3]]), Rationals());
        Gs +:= ee((a*(nv*Qr, nv)/2 + (vs[i]*Qr, nv)/dn)/c);
    end for;
    vf := ee(1/8) * CC!c^(-3/2) / (D*N*Sqrt(CC!2)) * ee(d*Qe[i]/c) * Gs;
    printf "eta %o: Q = %o  v_fft = %o  v_grp = %o  v_formula = %o\n", i, Qe[i], rnd(v2[i]), rnd(vg), rnd(vf);
    if Abs(vg) gt 10^-10 then Include(~r12, rnd(v2[i]/vg)); Include(~r32, rnd(vf/vg)); end if;
end for;
printf "ratios fft/grp: %o\nratios formula/grp: %o\n", r12, r32;
quit;
