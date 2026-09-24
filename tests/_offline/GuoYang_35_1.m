import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^35(1)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-7, Infinity()>,
<-8, 1>,
<-15, -3>,
<-28, 0>,
<-35, -7>,
<-43, 9>,
<-60, -1/3>,
<-67, 1/9>,
<-163, 81/25>,
<-235, -47/9>,
<-280, -9/7>,
<-315, -1/7>,
<-427, -175/9>,
<-595, 17/9>,
<-1435, -369/7>
];
test_gy_table(35, 1, gy);
