import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^39(1)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-7, Infinity()>,
<-15, -1>,
<-19, -3>,
<-24, 1>,
<-28, 0>,
<-60, -5>,
<-67, 7/3>,
<-91, -1/3>,
<-123, 1/5>,
<-163, -57/35>,
<-195, 5>,
<-267, -47/5>,
<-312, -17/7>,
<-403, -75/17>
];
test_gy_table(39, 1, gy);
