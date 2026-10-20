import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^22(3)", primary hauptmodule column. See GuoYangCheck.m for the method.
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
test_gy_table(22, 3, gy : Force := [-3, -132, -168, -267, -312, -372, -408, -627, -660, -708]);   // forced: the published points the search never offers -- conductor points, and those a level prime shares with d_0 (the ramified-level values; none is excluded)
