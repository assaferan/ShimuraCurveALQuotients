import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^74(1)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-8, Infinity()>,
<-19, -2>,
<-20, -1>,
<-24, -3>,
<-43, 0>,
<-52, 1>,
<-88, -5>,
<-148, -9/4>,
<-163, 10>
];
test_gy_table(74, 1, gy);
