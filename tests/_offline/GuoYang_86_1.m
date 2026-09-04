import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^86(1)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-4, Infinity()>,
<-11, 2>,
<-24, 3>,
<-40, 1>,
<-43, 9/4>,
<-52, 5/2>,
<-67, 0>,
<-232, -5>
];
test_gy_table(86, 1, gy);
