import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^58(1)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-3, 27/19>,
<-8, Infinity()>,
<-11, 11/19>,
<-19, 1>,
<-27, 3/19>,
<-40, -5/19>,
<-43, 43/19>,
<-148, 25/19>,
<-163, 163/475>,
<-232, 0>
];
test_gy_table(58, 1, gy);
