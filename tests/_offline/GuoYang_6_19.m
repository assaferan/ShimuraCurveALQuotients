import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^6(19)", primary hauptmodule column. See GuoYangCheck.m for the method.
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
test_gy_table(6, 19, gy : Force := [-19, -228]);   // forced: the published points the search never offers -- conductor points, and those a level prime shares with d_0 (the ramified-level values; none is excluded)
