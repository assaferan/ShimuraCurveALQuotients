import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^6(17)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-4, 0>,
<-19, 1>,
<-43, 9>,
<-51, Infinity()>,
<-52, -1>,
<-67, 1/4>,
<-84, -3>,
<-120, -1/3>,
<-123, 1/9>,
<-132, 3>,
<-408, -16/3>
];
test_gy_table(6, 17, gy);
