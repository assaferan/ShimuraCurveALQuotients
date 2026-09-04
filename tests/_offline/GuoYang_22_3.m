import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^22(3)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-3, 1>,
<-11, Infinity()>,
<-20, 5/4>,
<-132, 0>,
<-168, 27/28>,
<-267, 169/196>,
<-312, 25/52>,
<-372, 31/4>,
<-408, 18/17>,
<-627, -11/16>,
<-660, 45/44>,
<-708, 675/676>
];
test_gy_table(22, 3, gy);
