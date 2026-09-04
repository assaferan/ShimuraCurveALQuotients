import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^6(31)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-3, Infinity()>,
<-24, 0>,
<-43, 4/3>,
<-52, 16/3>,
<-84, 16/9>,
<-88, 16/27>,
<-120, 8/3>,
<-123, -16/9>,
<-148, 64/27>,
<-168, 8/9>,
<-228, -16/3>,
<-232, 1/3>,
<-372, 1>,
<-403, 52/27>
];
test_gy_table(6, 31, gy);
