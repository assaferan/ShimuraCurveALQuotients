import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^146(1)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-11, Infinity()>,
<-20, 1>,
<-40, -1>,
<-43, 0>,
<-52, 3>,
<-88, 2>,
<-232, 5>
];
test_gy_table(146, 1, gy);
