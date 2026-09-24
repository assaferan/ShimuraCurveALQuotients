import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^14(3)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-8, 0>,
<-11, 1>,
<-35, -1/7>,
<-51, 1/9>,
<-84, Infinity()>,
<-120, -1/27>,
<-123, 25/9>,
<-168, -2/9>,
<-228, -25/27>,
<-267, 25/1521>,
<-312, 49/117>
];
test_gy_table(14, 3, gy);
