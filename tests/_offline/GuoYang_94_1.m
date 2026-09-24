import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^94(1)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-3, -1>,
<-4, -9>,
<-8, Infinity()>,
<-24, -5>,
<-27, -7>,
<-148, -17>,
<-235, -13/5>
];
test_gy_table(94, 1, gy);
