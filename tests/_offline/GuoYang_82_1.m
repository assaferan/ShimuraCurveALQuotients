import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^82(1)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-3, 1/2>,
<-11, 0>,
<-19, Infinity()>,
<-24, 1>,
<-27, -1>,
<-52, 1/3>,
<-67, 2/3>,
<-88, -1/2>,
<-123, 1/4>,
<-232, 7/13>
];
test_gy_table(82, 1, gy);
