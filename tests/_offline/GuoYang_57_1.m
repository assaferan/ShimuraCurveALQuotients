import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^57(1)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-4, Infinity()>,
<-7, -1>,
<-16, 0>,
<-19, -1/3>,
<-24, 1>,
<-28, 1/3>,
<-43, -3>,
<-123, -1/2>,
<-163, -11/6>,
<-267, -11>
];
test_gy_table(57, 1, gy);
