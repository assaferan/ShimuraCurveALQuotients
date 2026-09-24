import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^87(1)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-3, Infinity()>,
<-12, 0>,
<-15, -3>,
<-19, 1>,
<-43, 9>,
<-48, -1>,
<-60, -1/3>,
<-147, -9/7>,
<-435, -5/3>
];
test_gy_table(87, 1, gy);
