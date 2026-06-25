# ShimuraCurveALQuotients
 Code for Atkin-Lehner quotients of Shimura Curves

# Loading the data

 In order to load the data in Magma, use the following Magma code.
 ```
 > AttachSpec("ShimuraQuotients.spec");
 > curves := GetHyperellipticCandidates();
 ```
 Each curve `X` stores it discriminant, level, and Atkin-Lehner subgroup,
 as attributes `D, N, W`, as well as its genus `g`.
 In addition, it stores the attribute `IsSubhyp` which, when assigned, 
 determines whether the curve is subhyperelliptic or not. 
 For example, to get all the curves that we have determined to be
 subhyperelliptic, use
 ```
 > subhyp := [X : X in curves | assigned X`IsSubhyp and X`IsSubhyp];
 > #subhyp;
 5846
 > subhyp_not_modular := [X : X in curves | assigned X`IsSubhyp and X`IsSubhyp and X`D ne 1];
 > #subhyp_not_modular;
 5367
 ```
 For those that have been determined as non-hyperelliptic, use
 ```
 > nonhyp := [X : X in curves | assigned X`IsSubhyp and not X`IsSubhyp];
 > #nonhyp;
 11706
 > nonhyp_not_modular := [X : X in curves | assigned X`IsSubhyp and not X`IsSubhyp and X`D ne 1];
 > #nonhyp_not_modular;
 10754
 ```
 For the curves that we were unable to determine, use
 ```
 > unknown := [X : X in curves | not assigned X`IsSubhyp];
 > #unknown;
 827
 > unknown_not_modular := [X : X in curves | not assigned X`IsSubhyp and X`D ne 1];
 > #unknown_not_modular;
 784
 ```
# The sets in Theorem 1.2.2
 The sets $$\mathcal{S}_h$$ and $$\mathcal{S}_u$$ from the paper can be obtained as follows.
 ```
 > S_h := [X : X in curves | assigned X`IsSubhyp and X`IsHyp];
 > #S_h;
 3423;
 > S_u := S_h cat unknown_not_modular;
 > #S_u;
 4207
 ```
# Additional functionality
 The package also includes some functionality to assist with more complicated queries,
 such as the intrinsic `IsStarCurve`, and each curve has attributes `Covers` and `CoveredBy`,
 which list the id numbers `CurveID` of the corresponding curves.

 Each curve has an attribute `TestInWhichProved`. When the curve has been classified 
 as either sub-hyperelliptic or non-hyperelliptic, the test proving it (and the relevant parameters)
 are stored there. For example
 ```
 > X := subhyp_not_modular[1]; X;
 Shimura quotient of level 1, and discriminant 6 by Atkin-Lehners { 1 }
 > X`TestInWhichProved;
 HyperellipticALInvolution to curve #1478
 > X := nonhyp_not_modular[1]; X;
 Shimura quotient of level 23, and discriminant 6 by Atkin-Lehners { 1 }
 > X`TestInWhichProved;
 ALFixedPointsOnQuotient, W_46 has 8 fixed points
 ```
 The package also produces models of curves via the method of Guo--Yang.

 For example, the models in Tables 3 and 4 can be reproduced by running the verification files VERIFICATION.m and VERIFICATION_genus0.m.
 More generally, to produce a desired model above a genus 0 star quotient one can run the following code.
 ```
 assert exists(Xstar){X : X in curves | X`D eq 26 and X`N eq 1 and IsStarCurve(X)};
 covers, ws := AllEquationsAboveCovers(Xstar, curves);
```

 This returns the equations for the fovers of X0*(26,1) in an asssociative array whose keys are the CurveIDs of curves.
