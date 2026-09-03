import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^62(1)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-4, Infinity()>,
<-8, 0>,
<-19, 1>,
<-20, -1>,
<-40, -1/2>,
<-67, 1/4>,
<-163, 25>,
<-372, -4/9>,
<-403, 13/4>
];
test_gy_table(62, 1, gy);
