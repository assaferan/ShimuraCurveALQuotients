import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^6(37)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-3, 0>,
<-4, 1>,
<-40, 9>,
<-67, 9/25>,
<-84, -9/7>,
<-120, -27/5>,
<-123, -3>,
<-132, 27/11>,
<-148, Infinity()>,
<-232, 81/49>,
<-312, -3/13>,
<-408, 9/17>,
<-555, -27/37>
];
test_gy_table(6, 37, gy);
