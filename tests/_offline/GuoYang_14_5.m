import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^14(5)", primary hauptmodule column. See GuoYangCheck.m for the method.
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
test_gy_table(14, 5, gy : Force := [-35, -120, -235, -280, -340, -420, -520, -840]);   // forced: the published points the search never offers -- conductor points, and those a level prime shares with d_0 (the ramified-level values; none is excluded)
