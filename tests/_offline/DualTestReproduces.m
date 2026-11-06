// The intrinsics of DualTest.m reproduce, from the saved bases, every certified relation and every
// refused control of the dual tests of tests/_offline/Hyperelliptic*.m (15_38, 21_58, 77_10,
// 39_34, 3390_1, 6_649, 6_385, 714_5, 4130_1), and the dimensions of the test spaces those tests
// record.  The bases are read from DualTestBasisDirectory() (DUALW_CACHE, default ~/dualw); the
// four below level 3000 are computed there if missing (minutes each), the others take days and
// must be present (tools/proof-4130-1/basisS.m).  CASES=label,... in the environment runs some cases only.  Each case: <label, D, N, M, B, variant,
// eigenspace <p, a> or <>, dimension of the space or 0, certified targets, refused targets>.

dif := func< X, Y | X cat [ <t[1], -t[2]> : t in Y ] >;
F2 := [<-16, 1>, <-120, 1>]; F4 := [<-88, 1>, <-100, 1>]; F5 := [<-112, 1>, <-372, 1>]; F6 := [<-168, 1>, <-228, 1>];
G2 := [<-16, 1>, <-100, 1>]; G3 := [<-280, 1>, <-660, 1>]; G7 := [<-760, 1>, <-1540, 1>]; G8 := [<-840, 1>, <-1320, 1>]; GQ := [<-4180, 1>];
H := [<-84, 1>, <-132, 1>]; H1 := [<-15, 1>, <-52, 1>]; H7 := [<-1092, 1>, <-372, 1>]; HQ1 := [<-3432, 1>]; HQ2 := [<-1560, 1>];
K := [ [<-35, 1>, <-532, 1>], [<-163, 1>, <-595, 1>], [<-315, 1>, <-520, 1>], [<-340, 1>, <-420, 1>] ];
cases := [* 
  <"15_38", 15, 38, 1140, 4000, "sharp", <7, -4>, 1,
     [ [<-228, 1>, <-760, 1>, <-12, -1>, <-40, -1>] ],
     [ [<-228, 1>, <-12, 1>, <-760, -1>, <-40, -1>], [<-228, 1>, <-40, 1>, <-760, -1>, <-12, -1>],
       [<-228, 1>, <-760, -1>], [<-12, 1>, <-40, -1>], [<-532, 1>, <-40, -1>] ]>,
  <"21_58", 21, 58, 2436, 4000, "plus2", <>, 3,
     [ dif(F2, F4), dif(F2, F5), dif(F2, F6), dif(F4, F5) ],
     [ [<-16, 1>, <-120, -1>], [<-88, 1>, <-100, -1>], [<-16, 1>, <-88, -1>], [<-112, 1>, <-372, -1>], [<-16, 1>, <-168, -1>] ]>,
  <"77_10", 77, 10, 1540, 8000, "plus2", <>, 3,
     [ dif(G2, G3), dif(G2, G7), dif(G2, G8), dif(G2, GQ) ],
     [ [<-16, 1>, <-100, -1>], [<-280, 1>, <-660, -1>], [<-16, 1>, <-280, -1>], [<-760, 1>, <-1540, -1>], [<-840, 1>, <-1320, -1>] ]>,
  <"39_34", 39, 34, 2652, 4000, "plus2", <>, 3,
     [ dif(H7, H), dif(HQ1, H), dif(HQ2, H) ],
     [ [<-84, 1>, <-132, -1>], [<-1092, 1>, <-372, -1>], [<-84, 1>, <-372, -1>] ]>,
  <"39_34", 39, 34, 2652, 4000, "odd", <>, 3,
     [ dif(H1, H), dif(H7, H) ],
     [ [<-15, 1>, <-52, -1>], [<-15, 1>, <-84, -1>] ]>,
  <"3390_1", 3390, 1, 6780, 1100, "sharp", <>, 51,
     [ [<-3, 1>, <-40, 1>, <-43, -1>, <-67, -1>], [<-3, 1>, <-40, 1>, <-123, -1>, <-168, -1>], [<-43, 1>, <-67, 1>, <-123, -1>, <-168, -1>] ],
     [ [<-3, 1>, <-40, -1>], [<-3, 1>, <-43, 1>, <-40, -1>, <-67, -1>] ]>,
  <"6_649", 6, 649, 7788, 1200, "sharp", <>, 53,
     [ [<-24, 1>, <-132, 1>, <-40, -1>, <-187, -1>], [<-24, 1>, <-132, 1>, <-52, -1>, <-660, -1>], [<-40, 1>, <-187, 1>, <-52, -1>, <-660, -1>] ],
     [ [<-24, 1>, <-132, -1>], [<-24, 1>, <-40, 1>, <-132, -1>, <-187, -1>] ]>,
  <"6_385", 6, 385, 4620, 4000, "sharp", <>, 1,
     [ [<-40, 1>, <-1540, 1>, <-9240, -1>] ], [ [<-40, 1>, <-1540, -1>], [<-40, 2>, <-9240, -1>] ]>,
  <"714_5", 714, 5, 7140, 4000, "sharp", <>, 1,
     [ [<-340, 1>, <-840, 1>, <-14280, -1>] ], [ [<-340, 1>, <-840, -1>], [<-340, 2>, <-14280, -1>] ]>,
  <"4130_1", 4130, 1, 8260, 4000, "sharp", <>, 3,
     [ dif(K[i], K[j]) : i, j in [1..4] | i lt j ],
     [ [<-35, 1>, <-532, -1>], [<-35, 1>, <-163, 1>, <-532, -1>, <-595, -1>] ]>
*];
only := GetEnv("CASES") ne "" select Split(GetEnv("CASES"), ",") else [];
for c in cases do
    label, D, N, M, B, variant, eig, dim, yes, no := Explode(c);
    if #only gt 0 and label notin only then continue; end if;
    printf "%o (%o, level %o, precision %o)...", label, variant, M, B;
    t0 := Realtime();
    A0 := ReadWeightThreeHalvesBasis(M, B : Compute := M lt 3000);
    T := DualTestSpace(D, N, A0 : Variant := variant);
    if #eig gt 0 then T := RestrictDualTestSpace(T, eig[1], eig[2]); end if;
    assert #T`Space eq dim;
    assert forall{t : t in yes | CertifiesTorsion(T, t)};
    assert forall{t : t in no | not CertifiesTorsion(T, t)};
    printf " dim %o, %o certified, %o refused (%o s)\n", dim, #yes, #no, Realtime() - t0;
end for;
