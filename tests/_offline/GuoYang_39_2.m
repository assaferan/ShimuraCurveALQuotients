import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^39(2)", primary hauptmodule column. See GuoYangCheck.m for the method.
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
test_gy_table(39, 2, gy : Force := [-24, -28, -52, -60, -84, -132, -148, -228, -232, -312, -372, -408, -520, -708, -1092]);   // forced: the published points the search never offers -- conductor points, and those a level prime shares with d_0 (the ramified-level values; none is excluded)
