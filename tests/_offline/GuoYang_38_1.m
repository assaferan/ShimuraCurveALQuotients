import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^38(1)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-4, Infinity()>,
<-11, 1>,
<-19, 0>,
<-20, -1>,
<-24, -3>,
<-43, 9>,
<-163, 4/81>,
<-228, -19/9>,
<-232, 29/81>,
<-532, -9/49>,
<-760, -171/16>
];
test_gy_table(38, 1, gy);
