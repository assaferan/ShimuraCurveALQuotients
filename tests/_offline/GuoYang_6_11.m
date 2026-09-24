import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^6(11)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-19, -1>,
<-24, 0>,
<-40, -1/9>,
<-43, -1/49>,
<-51, -1/17>,
<-52, -1/25>,
<-84, 1/7>,
<-88, Infinity()>,
<-120, 3/5>,
<-123, -9/41>,
<-132, 1>
];
test_gy_table(6, 11, gy);
