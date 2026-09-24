import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^22(5)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-4, 1>,
<-11, Infinity()>,
<-20, 0>,
<-115, 5/4>,
<-235, 5>,
<-280, -1/7>,
<-520, 5/13>,
<-660, -5/4>,
<-715, -5/11>,
<-760, -5/76>
];
test_gy_table(22, 5, gy);
