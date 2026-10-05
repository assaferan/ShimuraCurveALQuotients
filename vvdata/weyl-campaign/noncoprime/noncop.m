AttachSpec("ShimuraQuotients.spec");
_ := ClassNumberLU(-4);
import "tests/_offline/GuoYangCheck.m" : test_gy_table;
SetColumns(0);
printf "\n=== X0^10(11): non-coprime table points [-88, -132, -187, -660, -715], of which level prime divides d_0: [-88, -132, -187, -660, -715]\n";
gy := [
<-8, Infinity()>,
<-35, 7/8>,
<-40, 1>,
<-43, -1/8>,
<-52, 9/8>,
<-88, 0>,
<-120, 3/8>,
<-132, 11/8>,
<-187, -9/8>,
<-340, 17/8>,
<-660, 33/32>,
<-715, 99/104>
];
try test_gy_table(10, 11, gy : Force := [-88, -132, -187, -660, -715]); printf "\n"; catch e printf "\n  FAILED: %o\n", e`Object; end try;
printf "\n=== X0^10(13): non-coprime table points [-52, -195, -312], of which level prime divides d_0: [-52, -195, -312]\n";
gy := [
<-3, 1>,
<-35, 5>,
<-40, Infinity()>,
<-43, 9>,
<-52, 0>,
<-88, 25/9>,
<-120, -15>,
<-195, -5/3>,
<-235, 5/9>,
<-312, -13/3>
];
try test_gy_table(10, 13, gy : Force := [-52, -195, -312]); printf "\n"; catch e printf "\n  FAILED: %o\n", e`Object; end try;
printf "\n=== X0^10(23): non-coprime table points [-115], of which level prime divides d_0: [-115]\n";
gy := [
<-20, 0>,
<-40, Infinity()>,
<-43, 5/4>,
<-67, 5>,
<-88, 5/16>,
<-115, 1>,
<-120, 5/8>,
<-148, 5/9>,
<-235, 25/16>,
<-520, -5/8>
];
try test_gy_table(10, 23, gy : Force := [-115]); printf "\n"; catch e printf "\n  FAILED: %o\n", e`Object; end try;
printf "\n=== X0^14(3): non-coprime table points [-51, -84, -120, -123, -168, -228, -267, -312], of which level prime divides d_0: [-51, -84, -120, -123, -168, -228, -267, -312]\n";
gy := [
<-8, 0>,
<-11, 1>,
<-35, -1/7>,
<-51, 1/9>,
<-84, Infinity()>,
<-120, -1/27>,
<-123, 25/9>,
<-168, -2/9>,
<-228, -25/27>,
<-267, 25/1521>,
<-312, 49/117>
];
try test_gy_table(14, 3, gy : Force := [-51, -84, -120, -123, -168, -228, -267, -312]); printf "\n"; catch e printf "\n  FAILED: %o\n", e`Object; end try;
printf "\n=== X0^14(5): non-coprime table points [-35, -120, -235, -280, -340, -420, -520, -840], of which level prime divides d_0: [-35, -120, -235, -280, -340, -420, -520, -840]\n";
gy := [
<-4, Infinity()>,
<-11, 1>,
<-35, 0>,
<-84, -1>,
<-91, -7>,
<-120, 5>,
<-235, 25/81>,
<-280, 5/16>,
<-340, -25>,
<-420, 5/9>,
<-520, 5/81>,
<-840, -35/9>
];
try test_gy_table(14, 5, gy : Force := [-35, -120, -235, -280, -340, -420, -520, -840]); printf "\n"; catch e printf "\n  FAILED: %o\n", e`Object; end try;
printf "\n=== X0^15(2): non-coprime table points [-12, -28, -40, -48, -52, -60, -88, -120, -132, -148, -168, -228, -232, -240, -280, -312, -340, -372, -408, -420, -520, -660, -708, -760, -840], of which level prime divides d_0: [-40, -52, -88, -120, -132, -148, -168, -228, -232, -280, -312, -340, -372, -408, -420, -520, -660, -708, -760, -840]\n";
gy := [
<-7, 1/4>,
<-12, Infinity()>,
<-15, 5/4>,
<-28, 9/4>,
<-40, 0>,
<-48, -1/4>,
<-52, 1>,
<-60, -1/12>,
<-88, 4>,
<-120, 2>,
<-132, -1>,
<-148, 1/25>,
<-168, 2/3>,
<-228, -1/9>,
<-232, 144/121>,
<-240, -25/12>,
<-280, 10>,
<-312, 2/25>,
<-340, 9/17>,
<-372, -31/9>,
<-408, 68/25>,
<-420, 5/3>,
<-520, -8/121>,
<-660, -5/11>,
<-708, -841/121>,
<-760, 450/529>,
<-840, 40/27>
];
try test_gy_table(15, 2, gy : Force := [-12, -28, -40, -48, -52, -60, -88, -120, -132, -148, -168, -228, -232, -240, -280, -312, -340, -372, -408, -420, -520, -660, -708, -760, -840]); printf "\n"; catch e printf "\n  FAILED: %o\n", e`Object; end try;
printf "\n=== X0^21(2): non-coprime table points [-4, -16, -28, -60, -84, -100, -112, -120, -148, -168, -228, -232, -280, -312, -372, -408, -420, -532, -708, -840], of which level prime divides d_0: [-4, -16, -84, -100, -120, -148, -168, -228, -232, -280, -312, -372, -408, -420, -532, -708, -840]\n";
gy := [
<-4, Infinity()>,
<-7, -7>,
<-15, -5/3>,
<-16, 1>,
<-28, 1/9>,
<-60, 9>,
<-84, -3>,
<-100, 1/5>,
<-112, 25>,
<-120, -1/3>,
<-148, 37/9>,
<-168, 0>,
<-228, -25/3>,
<-232, -32>,
<-280, -35/9>,
<-312, -8/3>,
<-372, -3/4>,
<-408, -75>,
<-420, 21>,
<-532, -19/4>,
<-708, 25/48>,
<-840, -16/3>
];
try test_gy_table(21, 2, gy : Force := [-4, -16, -28, -60, -84, -100, -112, -120, -148, -168, -228, -232, -280, -312, -372, -408, -420, -532, -708, -840]); printf "\n"; catch e printf "\n  FAILED: %o\n", e`Object; end try;
printf "\n=== X0^22(3): non-coprime table points [-3, -132, -168, -267, -312, -372, -408, -627, -660, -708], of which level prime divides d_0: [-3, -132, -168, -267, -312, -372, -408, -627, -660, -708]\n";
gy := [
<-3, 1>,
<-11, Infinity()>,
<-20, 5/4>,
<-132, 0>,
<-168, 27/28>,
<-267, 169/196>,
<-312, 25/52>,
<-372, 31/4>,
<-408, 18/17>,
<-627, -11/16>,
<-660, 45/44>,
<-708, 675/676>
];
try test_gy_table(22, 3, gy : Force := [-3, -132, -168, -267, -312, -372, -408, -627, -660, -708]); printf "\n"; catch e printf "\n  FAILED: %o\n", e`Object; end try;
printf "\n=== X0^22(5): non-coprime table points [-20, -115, -235, -280, -520, -660, -715, -760], of which level prime divides d_0: [-20, -115, -235, -280, -520, -660, -715, -760]\n";
gy := [
<-4, 1>,
<-11, Infinity()>,
<-20, 0>,
<-115, 5/4>,
<-235, 5>,
<-280, -1/7>,
<-520, 5/13>,
<-660, -5/4>,
<-715, -5/11>,
<-760, -5/76>
];
try test_gy_table(22, 5, gy : Force := [-20, -115, -235, -280, -520, -660, -715, -760]); printf "\n"; catch e printf "\n  FAILED: %o\n", e`Object; end try;
printf "\n=== X0^39(2): non-coprime table points [-24, -28, -52, -60, -84, -132, -148, -228, -232, -312, -372, -408, -520, -708, -1092], of which level prime divides d_0: [-24, -52, -84, -132, -148, -228, -232, -312, -372, -408, -520, -708, -1092]\n";
gy := [
<-7, -7>,
<-15, -15>,
<-24, Infinity()>,
<-28, 9>,
<-52, -9>,
<-60, 1>,
<-84, -3>,
<-132, -11>,
<-148, -1>,
<-228, -27>,
<-232, -9/2>,
<-312, 0>,
<-372, -25/3>,
<-408, -12>,
<-520, -10>,
<-708, -59>,
<-1092, -1/3>
];
try test_gy_table(39, 2, gy : Force := [-24, -28, -52, -60, -84, -132, -148, -228, -232, -312, -372, -408, -520, -708, -1092]); printf "\n"; catch e printf "\n  FAILED: %o\n", e`Object; end try;
printf "\n=== X0^6(11): non-coprime table points [-88, -132], of which level prime divides d_0: [-88, -132]\n";
gy := [
<-19, -1>,
<-24, 0>,
<-40, -1/9>,
<-43, -1/49>,
<-51, -1/17>,
<-52, -1/25>,
<-84, 1/7>,
<-88, Infinity()>,
<-120, 3/5>,
<-123, -9/41>,
<-132, 1>
];
try test_gy_table(6, 11, gy : Force := [-88, -132]); printf "\n"; catch e printf "\n  FAILED: %o\n", e`Object; end try;
printf "\n=== X0^6(17): non-coprime table points [-51, -408], of which level prime divides d_0: [-51, -408]\n";
gy := [
<-4, 0>,
<-19, 1>,
<-43, 9>,
<-51, Infinity()>,
<-52, -1>,
<-67, 1/4>,
<-84, -3>,
<-120, -1/3>,
<-123, 1/9>,
<-132, 3>,
<-408, -16/3>
];
try test_gy_table(6, 17, gy : Force := [-51, -408]); printf "\n"; catch e printf "\n  FAILED: %o\n", e`Object; end try;
printf "\n=== X0^6(19): non-coprime table points [-19, -228], of which level prime divides d_0: [-19, -228]\n";
gy := [
<-3, 0>,
<-19, Infinity()>,
<-40, -1/4>,
<-51, 3/4>,
<-52, -9/4>,
<-67, -9/16>,
<-84, -3/4>,
<-88, -1>,
<-132, 1/4>,
<-148, -1/36>,
<-228, 1>
];
try test_gy_table(6, 19, gy : Force := [-19, -228]); printf "\n"; catch e printf "\n  FAILED: %o\n", e`Object; end try;
printf "\n=== X0^6(29): non-coprime table points [-232], of which level prime divides d_0: [-232]\n";
gy := [
<-4, 1>,
<-24, 0>,
<-51, 9/17>,
<-52, 9>,
<-67, 1/9>,
<-88, 9/25>,
<-120, -3/5>,
<-123, 9/41>,
<-132, 3/11>,
<-168, -9/7>,
<-228, 27/19>,
<-232, Infinity()>,
<-267, 81/89>
];
try test_gy_table(6, 29, gy : Force := [-232]); printf "\n"; catch e printf "\n  FAILED: %o\n", e`Object; end try;
printf "\n=== X0^6(31): non-coprime table points [-372, -403], of which level prime divides d_0: [-372, -403]\n";
gy := [
<-3, Infinity()>,
<-24, 0>,
<-43, 4/3>,
<-52, 16/3>,
<-84, 16/9>,
<-88, 16/27>,
<-120, 8/3>,
<-123, -16/9>,
<-148, 64/27>,
<-168, 8/9>,
<-228, -16/3>,
<-232, 1/3>,
<-372, 1>,
<-403, 52/27>
];
try test_gy_table(6, 31, gy : Force := [-372, -403]); printf "\n"; catch e printf "\n  FAILED: %o\n", e`Object; end try;
printf "\n=== X0^6(37): non-coprime table points [-148, -555], of which level prime divides d_0: [-148, -555]\n";
gy := [
<-3, 0>,
<-4, 1>,
<-40, 9>,
<-67, 9/25>,
<-84, -9/7>,
<-120, -27/5>,
<-123, -3>,
<-132, 27/11>,
<-148, Infinity()>,
<-232, 81/49>,
<-312, -3/13>,
<-408, 9/17>,
<-555, -27/37>
];
try test_gy_table(6, 37, gy : Force := [-148, -555]); printf "\n"; catch e printf "\n  FAILED: %o\n", e`Object; end try;
