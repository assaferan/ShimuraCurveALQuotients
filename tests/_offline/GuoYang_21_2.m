import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^21(2)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-4, Infinity()>,
<-7, -7>,
<-15, -5/3>,
<-16, 1>,
<-28, 1/9>,
<-60, 9>,
<-84, -3>,
<-100, 1/5>,
<-112, 25>,
<-120, -1/3>,
<-148, 37/9>,
<-168, 0>,
<-228, -25/3>,
<-232, -32>,
<-280, -35/9>,
<-312, -8/3>,
<-372, -3/4>,
<-408, -75>,
<-420, 21>,
<-532, -19/4>,
<-708, 25/48>,
<-840, -16/3>
];
test_gy_table(21, 2, gy);
