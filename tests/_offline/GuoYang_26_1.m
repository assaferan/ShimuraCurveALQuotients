import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^26(1)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-8, Infinity()>,
<-11, 1>,
<-19, 9>,
<-20, 5>,
<-24, -3>,
<-52, 0>,
<-67, 81/25>,
<-91, -9/7>,
<-148, 333/25>,
<-163, 4761/1225>,
<-232, -225/98>,
<-312, -75/2>,
<-403, -1521/775>,
<-520, 1521/245>
];
test_gy_table(26, 1, gy);
