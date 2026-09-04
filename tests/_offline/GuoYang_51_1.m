import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^51(1)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-3, Infinity()>,
<-7, 1>,
<-12, 0>,
<-24, -1/3>,
<-28, 1/9>,
<-51, -3>,
<-163, 1/16>,
<-187, -11/9>,
<-267, -25/3>,
<-408, -1/6>
];
test_gy_table(51, 1, gy);
