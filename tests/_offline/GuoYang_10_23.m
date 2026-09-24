import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^10(23)", primary hauptmodule column. See GuoYangCheck.m for the method.
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
test_gy_table(10, 23, gy);
