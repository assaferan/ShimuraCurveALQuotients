import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^6(29)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-4, 1>,
<-24, 0>,
<-51, 9/17>,
<-52, 9>,
<-67, 1/9>,
<-88, 9/25>,
<-120, -3/5>,
<-123, 9/41>,
<-132, 3/11>,
<-168, -9/7>,
<-228, 27/19>,
<-232, Infinity()>,
<-267, 81/89>
];
test_gy_table(6, 29, gy);
