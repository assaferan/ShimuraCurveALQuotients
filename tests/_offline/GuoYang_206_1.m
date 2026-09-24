import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^206(1)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-4, 1>,
<-8, Infinity()>,
<-19, 3>,
<-52, -1>,
<-163, 9>,
<-232, -3>
];
test_gy_table(206, 1, gy);
