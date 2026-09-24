import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^10(13)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-3, 1>,
<-35, 5>,
<-40, Infinity()>,
<-43, 9>,
<-52, 0>,
<-88, 25/9>,
<-120, -15>,
<-195, -5/3>,
<-235, 5/9>,
<-312, -13/3>
];
test_gy_table(10, 13, gy);
