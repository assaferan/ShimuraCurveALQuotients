import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^134(1)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-4, Infinity()>,
<-19, 2>,
<-24, 3>,
<-40, 1>,
<-67, 9/4>,
<-88, 5>,
<-148, 5/2>,
<-163, 0>
];
test_gy_table(134, 1, gy);
