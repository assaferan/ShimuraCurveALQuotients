import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^55(1)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-3, -3/2>,
<-12, Infinity()>,
<-15, 1>,
<-27, 0>,
<-60, -1>,
<-67, 8/3>,
<-88, -2>,
<-115, -4>,
<-163, -5/6>,
<-187, 1/2>,
<-235, -4/11>,
<-715, -4/3>
];
test_gy_table(55, 1, gy);
