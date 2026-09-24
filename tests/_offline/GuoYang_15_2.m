import "tests/_offline/GuoYangCheck.m" : test_gy_table;

// Guo-Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193), appendix table
// "CM-values of X_0^15(2)", primary hauptmodule column. See GuoYangCheck.m for the method.
gy := [
<-7, 1/4>,
<-12, Infinity()>,
<-15, 5/4>,
<-28, 9/4>,
<-40, 0>,
<-48, -1/4>,
<-52, 1>,
<-60, -1/12>,
<-88, 4>,
<-120, 2>,
<-132, -1>,
<-148, 1/25>,
<-168, 2/3>,
<-228, -1/9>,
<-232, 144/121>,
<-240, -25/12>,
<-280, 10>,
<-312, 2/25>,
<-340, 9/17>,
<-372, -31/9>,
<-408, 68/25>,
<-420, 5/3>,
<-520, -8/121>,
<-660, -5/11>,
<-708, -841/121>,
<-760, 450/529>,
<-840, 40/27>
];
test_gy_table(15, 2, gy);
