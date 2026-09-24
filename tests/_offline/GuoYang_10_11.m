import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^10(11)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-8, Infinity()>,
<-35, 7/8>,
<-40, 1>,
<-43, -1/8>,
<-52, 9/8>,
<-88, 0>,
<-120, 3/8>,
<-132, 11/8>,
<-187, -9/8>,
<-340, 17/8>,
<-660, 33/32>,
<-715, 99/104>
];
test_gy_table(10, 11, gy);
