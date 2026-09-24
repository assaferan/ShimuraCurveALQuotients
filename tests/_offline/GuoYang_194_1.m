import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^194(1)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-19, Infinity()>,
<-20, 1>,
<-40, -1>,
<-52, -3>,
<-67, 0>,
<-148, -1/3>,
<-232, 5>
];
test_gy_table(194, 1, gy);
