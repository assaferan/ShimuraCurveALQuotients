a:=-13; b:=6; Eichler:=1;
NMax:=70;

A<i,j,k>:=QuaternionAlgebra<Rationals()|a,b>;
D:=Discriminant(A);
L:=AssociativeArray();
DualL:=AssociativeArray();

if (D eq 6) and (Eichler eq 11) then
  divs := [[1,3,4,5,6,9,11,12,22], [2,5,6,9,22,33], [2,3,4,5,6,7,11,12,33,66]];
  wts := [[-2,-2,-2,-2,2,-2,-2,2,-2], [-6,2,-4,2,2,2], [-10,-2,-2,-2,-6,6,-2,2,2,2]];
  divs1 := [[],[],[]];
  wts1 := [[],[],[]];
  offs := [[],[],[]];
  offsets := [[],[],[]];
end if;

/* D=119 (-7,17) */
if D eq 119 then
  divs:=[
    [ 1, 2, 3, 5, 6, 7, 8, 9, 10, 12, 13, 15, 16, 17, 18, 19, 20, 21, 24, 25, 26, 27, 30, 31, 32, 33, 34, 36, 38, 40, 41, 42, 43, 45, 47, 49, 147 ],
    [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 15, 16, 17, 18, 19, 20, 21, 24, 25, 26, 27, 28, 30, 31, 32, 33, 34, 35, 38, 40, 41, 42, 43, 45, 47, 49, 50, 52, 147 ],
    [ 1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 13, 15, 16, 17, 18, 19, 20, 21, 24, 25, 26, 27, 30, 31, 32, 33, 34, 35, 36, 38, 40, 41, 42, 43, 45, 47, 49, 50, 52, 119, 147]];
  wts:=[
    [ 1, 5/2, -4, 1, 1/2, 2, 3, -3, 7/2, -1, -9/4, 7/2, 5/2, 3/2, -2, 2, 1/2, -9/2, -2, 2, -7/4, 3/2, -1/2, -1, 1, -3/4, 1/2, 1/2, -1, 5/2, 3/2, 1/2, -1/2, 1, -3/4, -1/2, 1 ],
    [ -1/8, 5/4, -2, 1/8, 3, -3/2, -1, 3/4, -5/8, 1, 1/2, -9/8, 5/4, 3/4, -1, -7/8, 7/8, 3, -3/4, -3, 1/2, -11/8, -2, 1, -3/8, -3/2, 3/8, -1/4, 1, -1/8, -5/8, -1/2, 1/2, 5/8, -1/2, 2, -1/2, -1/2, 1/8, -1/8, 3/2 ],
    [ 5/2, 25/2, -20, 1/2, 11, -9/2, 12, -27/2, 17/2, -3, -45/4, 35/2, 21/2, -11/2, -19/2, 11/2, 23/2, -37/2, -18, 8, -59/4, -11/2, -3, -11, 9/2, -13/4, 11/2, -1/2, 3/2, -11/2, -1/2, 5/2, 4, -7/2, 11, -17/4, -7/2, 1/2, -1/2, 2, 7 ]];
  divs1:=[[ 3/4, 7/4 ], [ 3/4, 7/4 ], [ 3/4, 7/4 ]];
  wts1:=[[ -2, -4 ], [ -2, -2 ], [ -14, -20 ]];
  offs:=[[7,17], [7,17], [7,17]];
  offsets:=[[-17/4,-11/4], [-5/2,-13/16], [-23/4,-25/2]];
end if;

if D eq 115 then
  divs:=[[ 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 14, 15, 17, 19, 20, 21, 22, 24, 25, 28, 29, 
30, 31, 38, 40, 43, 50, 75, 100 ],
  [ 1, 2, 4, 6, 9, 16, 24, 25, 26, 29, 31, 50, 75, 100 ]];
  wts:=[[ -3, 5/2, -1/2, -1/2, -1/2, -4, -1, -5, 1/2, -3/2, -1, -1/2, 3/2, 1/2, -3/2, 
-3/2, -1, -2, 1, 1/2, -2, 1, -3, 1/2, -1/2, 1/2, -5/2, -1/2, -1/2 ],
  [ -29/12, 5/2, -1/12, -13/3, -21/4, -13/12, -8/3, 5/6, 1, -11/4, -29/12, -13/6, 
-1, -1/2 ]];
  divs1:=[[ 3/4, 1, 3 ], [ 3/4, 1, 2, 3 ]];
  wts1:=[[ -4, -2, 2 ], [ -3, -7/3, 2/3, 1 ]];
  offs:=[[5,23],[2,5]]; offsets:=[[19/2,-1/2],[1,41/4]];
end if;

/* D=111 (-37,6) */
if D eq 111 then
  divs:=[
    [ 1, 2, 3, 4, 5, 8, 9, 10, 14, 15, 16, 18, 23, 25, 26, 27, 28, 29, 32, 33, 34, 35, 38, 40, 41, 44, 46, 48, 49, 80, 92 ],
    [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 17, 18, 21, 23, 25, 26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 38, 40, 41, 44, 45, 46, 47, 48, 49, 50, 72, 80, 89, 92 ],
    [ 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 14, 16, 17, 18, 20, 21, 23, 25, 26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 37, 38, 40, 41, 44, 45, 46, 47, 48, 49, 50, 72, 80, 89, 92 ],
    [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 17, 23, 25, 26, 27, 28, 29, 32, 33, 34, 35, 38, 40, 41, 44, 45, 46, 47, 48, 49, 54, 72, 89, 92, 111 ]];
  wts:=[
    [ -1/4, 5/4, -1/12, -1/24, -3/2, -1, 1/12, 1/12, -3/4, 1, -1/24, -1, 1/4, 1/12, 1/8, -1/6, -1/12, -3/4, -3/4, -1/3, 1/8, -1/2, -1/6, -1/8, 1/12, -1/24, 1/8, 1/6, 1/6, 1/2, -1 ],
    [ -1/4, 5/8, 1/12, -3/16, -5/4, 1, -1/12, -3/2, -1/24, -1/6, -5/24, -1/24, -3/8, 1/16, -1/3, -1/2, -3/8, -17/24, 1/12, 19/48, -1/8, -1/8, -3/8, -7/24, -17/24, -1/4, 3/16, 3/4, -1/12, 1/8, -1/48, -1/6, -7/48, 1/2, -7/48, -1/24, -1/8, -1/6, 1/6, 1/2, 1/4, -1/2, -1/2 ],
    [ 3/2, 5/4, 1/6, -3/8, -1/2, -1/6, -1, -1/12, 5/3, -5/12, -1/12, -3/4, 1/8, 4/3, 3, 2, 5/4, -41/12, 1/6, -29/24, -1/4, -1/4, -3/4, 17/12, -17/12, 3/2, 3/8, 3/2, -1/6, 2, 1/4, -1/24, -1/3, -7/24, 3, -7/24, -1/12, -1/4, -1/3, 1/3, 1, -3/2, -1, 1 ],
    [ -1, 5, 5/3, 11/6, -4, -2, 4, -6, 7/3, -5/3, 2, 4, -1, 11/6, 2, -3, -5/3, 5/2, -2/3, 5/3, -3, -3, 2/3, 1/2, 2, 4/3, -1/2, 1/3, -1/6, 4, 1/2, 2, 2/3, 2/3, 2, 2, -4, -2, 2 ]];
  divs1:=[[ 15/4 ], [ 3/2, 15/4 ], [ 3/2, 15/4 ], [ 15/4 ]];
  wts1:=[[ -2 ], [ 2, -1 ], [ 4, -2 ], [ -8 ]];
  offs:=[[3,37],[3,37],[3,37],[3,37]];
  offsets:=[[3,-1/12],[2/3,2/3],[-11/3,-5/3],[1,-19/3]];
end if;

/* D=95 (-5,38) */
if D eq 95 then
  divs:=[
    [ 1, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 21, 24, 25, 26, 27, 31, 36, 39, 40, 41, 44, 100, 125 ],
    [ 1, 3, 4, 5, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 21, 22, 25, 26, 27, 29, 31, 34, 36, 37, 39, 40, 41, 44, 100, 125 ],
    [ 1, 4, 6, 10, 11, 12, 13, 14, 15, 16, 18, 19, 21, 24, 25, 26, 27, 31, 36, 39, 40, 41, 44, 95 ],
    [ 1, 2, 3, 4, 5, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 19, 21, 22, 24, 25, 26, 27, 29, 31, 34, 36, 37, 39, 40, 41, 44, 100, 125, 190 ]];
  wts:=[
    [ 3/2, -2, -1/2, -2, 1, -1, 1/2, -3/2, -1/4, 1/4, -1, -1/4, -1, 1/4, -1, 1/4, -1, -1/2, -2, 1/4, 1/4, 1/2, -1/2, 1/4, 1/4, -1/2, -1/2, 1/2 ],
    [ -1/2, 1/2, 3, 5/2, 1/2, 1, 1/2, -1/2, -1/4, -1/4, 3, 1/4, -1, 1/4, 3, -1/4, -1/2, 1/2, 2, 1/4, 1/2, 3/4, 1/2, 1/2, -1/2, -1/2, 3/4, 3/4, 3/2, -1/2, -1/2 ],
    [ 6, -4, -4, 4, -2, -1, 1, 2, -1, -6, 3, 2, 3, -2, -2, -4, 3, 1, 2, -4, 1, 1, 2, 2 ],
    [ -5/4, -1, 1/4, 5/2, 1/4, 1/4, -1/2, 1/4, 3/4, -1/8, -1/8, 7/2, 1/8, 1/2, 1/8, 3/2, -9/8, -1/4, 1, 5/4, 1, 1/8, 1/4, 3/8, 1/4, 1/4, -1/4, 3/4, 3/8, 3/8, 3/4, -1/4, -1/4, 1 ]];
  divs1:=[[ 7/4 ], [ 5/4, 7/4 ], [ 7/4 ], [ 5/4, 7/4 ]];
  wts1:=[[ -2 ], [ 4, -2 ], [ -8 ], [ 2, -1 ]];
  offs:=[[5,19],[2,5,19],[2,5,19],[2,5,19]];
  offsets:=[[11/2,-3/8],[2,-3/2,-5/8],[8,6,-11/2],[1,-13/4,-5/16]];
end if;

/* D=93 (-31, 3) */
if D eq 93 then
  divs:=[
    [ 1, 2, 8, 9, 14, 20, 32, 35, 36, 41, 45, 71 ],
    [ 1, 2, 3, 4, 5, 6, 8, 9, 11, 13, 14, 17, 18, 20, 21, 23, 26, 27, 29, 30, 32, 34, 35, 36, 37, 38, 41, 42, 45, 71, 72 ],
    [ 1, 2, 4, 6, 7, 8, 9, 14, 15, 20, 22, 27, 30, 31, 32, 35, 36, 41, 45, 63, 71 ],
    [ 1, 2, 3, 5, 6, 8, 9, 11, 12, 15, 17, 18, 20, 21, 22, 23, 24, 26, 27, 29, 30, 34, 71, 93 ]];
  wts:=[
    [ -2, 2, 1/2, -1, -1/2, 2, 3/2, -1/2, 1, 1/2, 1, -1 ],
    [ -7/4, 19/24, 1/24, -3/4, 7/48, -1/6, 11/24, -1/4, 1/24, -1/8, 1/16, -1/12, 1/4, 31/24, -1/4, 1/24, 1/12, -1/24, 1/6, -1/24, 5/6, -1/8, -13/24, 3/4, -1/24, 7/48, 7/16, -5/24, 3/4, -3/4, 1/4 ],
    [ -10, 6, -2, 2, 2, 1, -2, 1, 2, 6, 2, 2, 2, 2, 3, -3, 2, 3, 6, -2, -2 ],
    [ -2, -4, -2, -2, 10, -2, -2, -2, 4, 4, -2, 2, -2, 6, 6, 2, 4, 2, 2, -4, 6, 4, -2, 2 ]];
  divs1:=[[], [7/4 ], [], []];
  wts1:=[[], [ 1 ], [], []];
  offs:=[[3],[3,31],[2,3,31],[3,31]];
  offsets:=[[-3], [-10/3,3/8],[6,-11,-4],[5,-17]];
end if;

/* D=87 (-3, 29) */
if D eq 87 then
  divs:=[
    [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 17, 18, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 38, 41, 44, 47, 50, 72, 90 ],
    [ 1, 2, 3, 5, 6, 8, 9, 10, 11, 12, 14, 16, 17, 18, 20, 22, 24, 25, 26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 38, 41, 44, 47, 50, 72, 90 ],
    [ 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 14, 16, 17, 18, 20, 22, 24, 25, 26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 38, 41, 44, 47, 50, 72, 87, 90 ]];
  wts:=[
    [ -3/4, -19/2, -2, 1/12, -5/4, -1/3, -1/6, 23, -5/12, 1, 1/12, 3, -31/2, -2, -17/12, 5/12, 1/12, -1/6, 2/3, -31/2, 4, -1/6, -3/2, -7/12, 7, -5/12, -7/12, -23/12, -5/12, -5/4, -6, -5, 10, 11, -3, -1 ],
    [ -1/2, -21/4, -2, -3/4, -1/4, 25/2, -1/4, 1/2, -1, 1, 1, 1/8, -29/4, -2, -3/4, 3/8, -1/8, 1/2, -31/4, 1, -1/8, -3/4, -3/8, 5/2, -1/4, -3/8, -1, -1/4, -5/8, -5/2, -2, 5, 6, -5/2, -1/2 ],
    [ -3, -47/2, -10, -1/2, 1/2, 2, 69, -3/2, 1, -4, 10, 3/4, -75/2, -8, -13/2, 17/4, -3/4, 3, -81/2, 10, -3/4, -9/2, -1/4, 17, -3/2, -9/4, -6, -3/2, -15/4, -17, -14, 28, 32, -11, 2, -1 ]];
  divs1:=[[ 3/4 ], [ 3/4 ], [ 3/4 ]];
  wts1:=[[ -8 ], [ -6 ], [ -36 ]];
  offs:=[[3,29], [3,29], [2,3,29]];
  offsets:=[[9/4,11/12], [33/8,7/16], [12,31/4,-11/8]];
end if;

/* D=69 (-3, 23) */
if D eq 69 then
/*  divs:=[
    [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15, 17, 18, 19, 20, 21, 22, 23, 26, 27, 28, 29, 32, 35, 36, 41, 47, 54 ],
    [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15, 17, 18, 19, 20, 21, 22, 23, 26, 27, 28, 29, 35, 36, 41, 47, 54 ],
    [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15, 17, 18, 19, 20, 21, 22, 23, 26, 27, 28, 29, 32, 35, 36, 41, 47, 54, 69 ]];
  wts:=[
    [ 5, 13, 4, -3, -7, 3, -5/6, -1/2, -5, 5/3, -47/6, -3, -5/2, 5/3, -4, 5/3, -5/2, 5/3, -5/6, -1, -5/2, -2, -5/6, -5/2, 1/2, -3/2, 3, 4, -9/2, -3 ],
    [ 9/2, 7, 1, -1/2, -2, 3/2, -1/3, -2, -5/2, 2/3, -11/6, -1, -1, 7/6, -3, 2/3, -1, 7/6, -5/6, -1, -3/2, -1, -1/3, -3/2, -1, 1/2, 2, -3/2, -3/2 ],
    [ 12, 25, 4, -2, -25/2, 6, -5/3, -6, -12, 4/3, -35/3, -9/2, -5, 23/6, -8, 10/3, -5, 23/6, -13/6, -2, -6, -4, -5/3, -3, 1, -4, 2, 9, -7, -6, 2 ]];
  divs1:=[[ 3/4 ], [ 3/4 ], [ 3/4 ]];
  wts1:=[[ -8 ], [ -6 ], [ -24 ]];
  offs:=[[3,23],[3,23],[2,3,23]];
  offsets:=[[33/2,-5/12],[11,-1/6],[8,40,-5/6]]; */
  divs:=[[ 1, 2, 3, 4, 5, 6, 8, 9, 11, 18, 23, 26, 27, 29, 35, 36, 41, 47, 54 ],
    [ 1, 2, 4, 5, 6, 8, 9, 11, 18, 23, 26, 29, 32, 36, 41, 54 ]];
  wts:=[[ 7, 14, 3, -2, 2, 1, -3, -3, 2, -5, -2, -2, -1, -5, -2, 2, 2, -3, -1 ],
    [ 7/2, 11/2, -1/2, 1, 3/2, -3/2, -3/2, 1, -3, -1, -2, -3/2, -1, 1/2, 1, -3/2 ]];
  divs1:=[[ 3/4, 1, 3/2 ], [ 3/4, 3/2, 9/4 ]];
  wts1:=[[ -8, -4, -12 ], [ -6, -6, -2 ]];
  offs:=[[3], [2,3]]; offsets:=[[17],[2,17/2]];
end if;

/* D=65 (-13,5) */
if D eq 65 then
  divs:=[[ 2, 6, 7, 11, 19, 21, 31, 50 ], [ 2, 6, 11, 19, 21, 31, 50 ],
    [ 1, 2, 4, 5, 6, 9, 10, 11, 13, 16, 17, 19, 21, 22, 31, 50 ]];
  wts:=[[ -1, 2, 1, -1, -2, -1, -1, 1 ],[ 16/12, -10/12, -4/12, 1/12, -4/12, 5/12,-4/12],
    [ 2, -14/3, -2, 2, 32/3, 2, -4, -16/3, 2, -2, -2, -29/3, -16/3, -2, -13/3, 14/3 ]];
  divs1:=[[ 7/4 ], [ 1/2, 5/4, 7/4 ],[ 1/2, 5/4, 7/4 ]];
  wts1:=[[ -2 ], [ -1, -3, -1],[ -12, -4, -4]];
  offs:=[[5],[2,5],[2,5,13]];
  offsets:=[[1],[1,11/12],[4,14/3,3]];
end if;

/* D=57 (-19,3) */
if D eq 57 then
  divs:=[
    [ 1, 5, 9, 11, 20, 36, 44 ],
    [ 1, 4, 5, 6, 9, 11, 17, 20, 23, 36, 44, 45 ],
    [ 1, 3, 4, 5, 9, 10, 11, 15, 17, 18, 19, 20, 36, 44, 45 ]];
  wts:=[
    [ -2, -1/6, -1, 5/6, -1/6, 1, 1 ],
    [ -5/3, -1, 31/12, 1, -1/3, 13/12, 5/2, -1/4, 1/2, 1, 7/6, -2/3 ],
    [ -22/3, 2, -4, 35/6, -2/3, 2, 35/6, 2, 4, 2, 2, -5/6, 4, 7/3, -4/3 ]];
  divs1:=[[], [], []];
  wts1:=[[], [], []];
  offs:=[[ 3 ], [ 2, 3 ], [ 2, 3, 19 ]];
  offsets:=[[ -11/3 ], [ 1, -37/6 ], [ 4, -50/3, -4]];
end if;

/* D=55 (-5,22) */
if D eq 55 then
  divs:=[
    [ 1, 4, 6, 7, 8, 9, 10, 11, 13, 14, 16, 25, 26, 75 ],
    [ 1, 3, 4, 6, 7, 8, 9, 10, 11, 13, 14, 16, 25, 26, 75 ],
    [ 1, 3, 4, 6, 7, 8, 9, 10, 11, 13, 14, 16, 25, 26, 55, 75 ]];
  wts:=[
    [ 8, -3, -7/2, 3/2, 3/2, 3, -3/2, -2, -3/2, -3, 1, -3, 4, -1 ],
    [ 21/4, -3/4, -1/4, -3/2, 1, 1, 2, -1, -1/2, -1, -5/2, 1/4, -7/4, 9/4, -1/4 ],
    [ 19, -1, -10, -8, 4, 2, 6, -2, -6, -4, -6, 2, -7, 10, 2, -3 ]];
  divs1:=[[ 3/4, 3], [ 3/4, 3, 15/4], [ 3/4, 3]];
  wts1:=[[ -4, 2], [ -3, 1, 1], [ -12, 4]];
  offs:=[[5,11],[2,11],[2,5,11]];
  offsets:=[[3,-3/4],[1,-1/2],[4,10,-2]];
  divs:=[
    [ 1, 4, 6, 7, 8, 9, 10, 11, 13, 14, 25, 26, 75 ],
    [ 1, 3, 4, 6, 7, 8, 9, 10, 11, 13, 14, 16, 22, 25, 26, 75 ],
    [ 1, 3, 4, 6, 7, 8, 9, 10, 11, 13, 14, 25, 26, 55, 75 ]];
  wts:=[
    [ 9, 3/2, 2, 1, 1, 9/2, -1, 3, -1, -9/2, -3, 3, 1 ],
    [ 15/4, 1/4, -5/4, -9/4, 3/4, 3/4, 5/4, -7/4, -1/2, -3/4, -5/4, 1/2, 1, -3/4, 3/2, -1/4 ],
    [ 21, -1, -1, 3, 3, 1, 9, -1, 4, -3, -9, -7, 8, 2, 1 ]];
  divs1:=[[ 3/4 ], [ 3/4 ], [ 3/4 ]];
  wts1:=[[ -4 ], [ -3 ], [ -12 ]];
  offs:=[[5,11],[5,11],[2,5,11]];
  offsets:=[[-9/2,-1/2],[1/4,-3/8],[4,-5,-3/2]];
end if;

/* D=51 (-3,17) */
if D eq 51 then
  divs:=[[ 3, 5, 8, 11, 14, 17, 20, 23, 27, 29, 45 ],
    [ 3, 5, 6, 8, 11, 14, 17, 20, 27, 29, 45 ],
    [ 2, 3, 4, 5, 8, 9, 11, 14, 17, 20, 23, 27, 29, 45, 51 ]];
  wts:=[[ -3, 10, -2, 3, 2, -2, 7, 2, 5, -2, -3 ],
    [ -3/2, 17/2, 1, -3/2, 1/2, 7/2, -3/2, 2, 3/2, 1/2, -5/2 ],
    [ -2, -7, -2, 28, -8, -2, 4, 10, -6, 9, 1, 7, -1, -9, 2 ]];
  divs1:=[[ 3/4 ], [ 3/4 ], [ 3/4 ]];
  wts1:=[[ -8 ], [ -6 ], [ -24 ]];
  offs:=[[3],[2,3],[2,3]];
  offsets:=[[-16],[2,-15/4-2],[8,-19]];
end if;

/* D=39 (-13,6) */
if D eq 39 then
  divs:=[[ 2, 5, 7, 8, 11, 18, 32 ], [ 2, 5, 6, 8, 11, 18, 32 ],
    [ 1, 2, 3, 8, 9, 10, 11, 14, 32, 39 ]];
  wts:=[[ -1/2, 1, 1, 3/2, 1/2, 1, -1 ],
    [ -1/4, 2/3, 1, -1/12, 1/4, -1/6, 1/3 ],
    [ 2, -1, 2, 1, -2, 2, 1, 2, -2, 2 ]];
  divs1:=[[ 7/4 ], [ 7/4 ], [ 7/4 ]];
  wts1:=[[ -2 ], [ -1 ], [ -4 ]];
  offs:=[[3],[2,3],[2,3,13]];
  offsets:=[[-2],[1,-1],[4,1,-4]];
end if;

/* D=35 (-7,5) */
if D eq 35 then
  divs:=[[ 1, 6, 7, 9, 10, 12, 25 ], [ 1, 2, 3, 5, 9, 10, 25 ],
    [ 1, 3, 4, 5, 6, 9, 10, 12, 25, 35 ]];
  wts:=[[ -6, -1, 2, 4, 2, 1, 2 ], [ -7/3, 1, -1/2, 1/2, 4/3, 1/2, 1/3 ],
    [ -38/3, -1, -2, 3, -1, 20/3, 3, -1, 8/3, 2 ]];
  divs1:=[[ 7/4 ], [ 7/4 ], [ 7/4 ]];
  wts1:=[[ -4 ], [ -2 ], [ -8 ]];
  offs:=[[5,7],[2,5,7],[2,5,7]];
  offsets:=[[-6,-2],[2,-4/3,1/4],[8,-26/3,-1/2]];
end if;

/* D=21 (-7,3) 1, (1+i)/2, j, (j+k)/2 */
if D eq 21 then
  divs:=[[ 1, 2, 3, 6, 7, 9 ], [ 1, 2, 3, 6, 9 ],
     [ 1, 2, 3, 5, 6, 9, 21 ]];
  wts:=[[ -2, 4, 4, 4, 2, 2 ], [ 2/3, 7/3, 1, 1, 4/3 ],
     [ -2/3, 2/3, 2, -2, 4, 2/3, 2 ]];
  divs1:=[[ 7/4 ], [ 7/4 ], [ 7/4 ]];
  wts1:=[[ -4 ], [ -2 ], [ -4 ]];
  offs:=[[3,7],[2,3,7],[2,3,7]];
  offsets:=[[-6,-5],[2,-23/6,-1/2],[4,-8/3,-4]];
end if;

/* D=15 (-3,5) */
if D eq 15 then
  divs:=[[2,3,8,2],[1,10,2],[1,2,8,15,2]];
  wts:=[[1,2,1,-4],[1,1,-3],[-2,1,-1,2,-6]];
  divs1:=[[3/4],[3/4],[3/4]];
  wts1:=[[-8],[-6],[-12]];
  offs:=[[],[2,3,5],[2,3,5]]; offsets:=[[],[-1,3,-3/2],[4,6,-1]];
end if;

/*
divs:=[]; wts:=[];
for i in [1..#divs2] do
  T:=Maximum(Maximum(divs0),Maximum(divs2[i]));
  LL:=AssociativeArray();
  for m in [1..T] do LL[m]:=0; end for;
  for j in [1..#divs0] do
    LL[divs0[j]]:=LL[divs0[j]]+wts0[j];
  end for;
  for j in [1..#divs2[i]] do
    LL[divs2[i][j]]:=LL[divs2[i][j]]+wts2[i][j];
  end for;
  M1:=[]; M2:=[];
  for j in [1..T] do
    if LL[j] ne 0 then
      Append(~M1,j); Append(~M2,LL[j]);
    end if;
  end for;
  Append(~divs,M1); Append(~wts,M2);
end for;
*/

P<z>:=PolynomialRing(Rationals());
LL<z>:=FieldOfFractions(P);
RR<X>:=PolynomialRing(LL);

Discs:={};
Optimals:=AssociativeArray();
ClassNos:=AssociativeArray();

QForm:=function(alpha,beta)
  return Trace(alpha*Conjugate(beta));
end function;

Convert:=function(alpha)
  return Vector(Rationals(),3,[Trace(alpha*i)/2/a,Trace(alpha*j)/2/b,
    -Trace(alpha*k)/2/a/b]);
end function;

OO:=QuaternionOrder(A,Eichler);
B:=Basis(OO);
for m in [2..4] do
  u:=Trace(B[m]);
  B[m]:=B[m]-Floor(u/2);
end for;
s:=0;

if Trace(B[2]) eq 1 then
  L[1]:=2*B[2]-1; s:=2;
else
  L[1]:=B[2];
end if;

if Trace(B[3]) eq 1 then
  if s eq 2 then
    L[2]:=B[3]-B[s];
  else
    L[2]:=2*B[3]-1; s:=3;
  end if;
else
  L[2]:=B[3];
end if;

if Trace(B[4]) eq 1 then
  if s ne 0 then
    L[3]:=B[4]-B[s];
  else
    L[3]:=2*B[4]-1;
  end if;
else
  L[3]:=B[4];
end if;

MM:=MatrixAlgebra(Rationals(),3)!0;

for m in [1..3] do
  LL:=Eltseq(L[m]);
  for n in [1..3] do
    MM[m,n]:=LL[n+1];
  end for;
end for;

LocalOptimal:=function(disc,p,n)
  x,y:=SquarefreeFactorization(disc);
  if x mod 4 eq 2 or x mod 4 eq 3 then
    x:=x*4; y:=y div 2;
  end if;
  s:=Valuation(y,p);
  if KroneckerSymbol(x,p) eq 1 then
    if n gt 2*s and s gt 0 then
      S:=2*(p+1)*p^(s-1);
    end if;
    if n gt 2*s and s eq 0 then
      S:=2;
    end if;
    if n eq 2*s and s gt 0 then
      S:=(p+2)*p^(s-1);
    end if;
    if n gt 0 and n lt 2*s and n mod 2 eq 0 then
      S:=(p+1)*p^((n div 2)-1);
    end if;
    if n lt 2*s and n mod 2 eq 1 then
      S:=2*p^((n-1) div 2);
    end if;
    if n eq 0 then
      S:=1;
    end if;
  end if;
  if KroneckerSymbol(x,p) eq -1 then
    if n gt 2*s then
      S:=0;
    end if;
    if n eq 2*s then
      S:=p^s;
    end if;
    if n gt 0 and n lt 2*s and n mod 2 eq 0 then
      S:=(p+1)*p^((n div 2)-1);
    end if;
    if n lt 2*s and n mod 2 eq 1 then
      S:=2*p^((n-1) div 2);
    end if;
    if n eq 0 then
      S:=1;
    end if;
  end if;
  if KroneckerSymbol(x,p) eq 0 then
    if n gt 2*s+1 then
      S:=0;
    end if;
    if n eq 2*s+1 then
      S:=p^s;
    end if;
    if n gt 0 and n lt 2*s+1 and n mod 2 eq 0 then
      S:=p^(n div 2)+p^((n div 2)-1);
    end if;
    if n gt 0 and n lt 2*s+1 and n mod 2 eq 1 then
      S:=2*p^((n-1) div 2);
    end if;
    if n eq 0 then
      S:=1;
    end if;
  end if;
  return S;
end function;

Proj:=function(alpha,S)
  beta:=alpha;
  for gamma in S do
    beta:=beta-QForm(beta,gamma)/QForm(gamma,gamma)*gamma;
  end for;
  return beta;
end function;

LCoefficients:=function(alpha)
  L:=Eltseq(alpha);
  W:=Solution(MM,Vector(Rationals(),3,[L[r]: r in [2..4]]));
  return W;
end function;

GramM:=MatrixAlgebra(Integers(),3)!0;
for m in [1..3] do
  for n in [1..3] do
    GramM[m,n]:=Trace(L[m]*Conjugate(L[n]));
  end for;
end for;

AA,U,V:=SmithForm(GramM);
for m in [1..3] do
  DualL[m]:=(U[m,1]*L[1]+U[m,2]*L[2]+U[m,3]*L[3])/AA[m,m];
end for;

XtoEvaluate:=function(alpha,m,eta)
  etaplus:=eta-Proj(eta,{alpha});
  u:=Rationals()!(etaplus/alpha);
  v:=Norm(alpha);
  return [(u+n)*alpha: n in [Ceiling(-Sqrt(m/v)-u)..Floor(Sqrt(m/v)-u)]]; 
end function;

SplitLattice:=function(Lplus)
  Mplus:=LCoefficients(Lplus);
  M:=Matrix(Rationals(),3,3,[LCoefficients(Proj(i,{Lplus})),
     LCoefficients(Proj(j,{Lplus})),LCoefficients(Proj(k,{Lplus}))]);
  B:=Basis(RowSpace(M));
  LL:=LatticeWithBasis(3,[1,0,0,0,1,0,0,0,1]);
  MM:=LatticeWithBasis(3,[B[1][1],B[1][2],B[1][3],B[2][1],
       B[2][2],B[2][3]]);
  B:=Basis(LL meet MM);
  R:=RSpace(Integers(),3);
  M:=Matrix(Integers(),3,3,[R!Mplus,R!B[1],R!B[2]]);
  return M;
end function;

GetDiag:=function(Sgcd,SS,p)
  if SS[1,1] mod p ne 0 then
    return Matrix(Rationals(),2,2,
      [SS[1,1]*Sgcd,0,0,Determinant(SS)/SS[1,1]*Sgcd]),
      Matrix(Rationals(),2,2,[1,0,-SS[1,2]/SS[1,1],1]);
  else
    if SS[2,2] mod p ne 0 then
      return Matrix(Rationals(),2,2,
        [Determinant(SS)/SS[2,2]*Sgcd,0,0,SS[2,2]*Sgcd]),
        Matrix(Rationals(),2,2,[1,-SS[1,2]/SS[2,2],0,1]);
    else
      if p ne 2 then
        A:=Matrix(Rationals(),2,2,[1,1,0,1]);
        T:=A*SS*Transpose(A);
        return Matrix(Rationals(),2,2,
          [T[1,1]*Sgcd,0,0,Determinant(T)/T[1,1]*Sgcd]),
          Matrix(Rationals(),2,2,[1,0,-T[1,2]/T[1,1],1])*A;
      else
        return SS*Sgcd, Matrix(Rationals(),2,2,[1,0,0,1]);
      end if;
    end if;
  end if;
end function;

Hmu:=function(p,mu)
  S:=[];
  for n in [1..2] do
    if Valuation(mu[n],p) ge 0 then
      Append(~S,n);
    end if;
  end for;
  return S;
end function;

Kmu:=function(p,L,mu)
  S:=10000000000000000;;
  if not (1 in Hmu(p,mu)) then
    S:=Valuation(L[1],p)+Valuation(mu[1],p);
  end if;
  if not (2 in Hmu(p,mu)) then
    u:=Valuation(L[2],p)+Valuation(mu[2],p);
    if u lt S then S:=u; end if;
  end if;
  return S;
end function;

Lmu:=function(p,L,mu,m)
  S:=[];
  for n in Hmu(p,mu) do
    ell:=Valuation(L[n],p);
    if ell lt m and IsOdd(ell-m) then
      Append(~S,n);
    end if;
  end for;
  return S;
end function;

lmu:=function(p,L,mu,m)
  return #Lmu(p,L,mu,m);
end function;

dmu:=function(p,L,mu,m)
  S:=m;
  for n in Hmu(p,mu) do
    ell:=Valuation(L[n],p);
    if ell lt m then
      S:=S+(ell-m)/2;
    end if;
  end for;
  return S;
end function;

vmu:=function(p,L,mu,m)
  if #Lmu(p,L,mu,m) lt 2 then
    S:=1;
  else
    S:=LegendreSymbol(-1,p);
  end if;
  for n in Lmu(p,L,mu,m) do
    eps:=L[n]/p^Valuation(L[n],p);
    S:=S*LegendreSymbol(Numerator(eps)*Denominator(eps),p);
  end for;
  return S;
end function;

fmu:=function(p,L,mu,u)
  c:=Valuation(u,p);
  e:=u/p^c;
  if #Lmu(p,L,mu,c+1) mod 2 eq 0 then
    return -1/z^2;
  else
    return LegendreSymbol(Numerator(e)*Denominator(e),p)/z;
  end if;
end function;

tmu:=function(p,L,mu,m)
  S:=m;
  for n in [1..2] do
    if not (n in Hmu(p,mu)) then
      S:=S-L[n]*mu[n]^2;
    end if;
  end for;
  return S;
end function;

Wmu:=function(disc,p,L,mu,m)
  S:=RR!0; T:=RR!0;
  c:=Valuation(tmu(p,L,mu,m),p);
  if Valuation(m-L[1]*mu[1]^2-L[2]*mu[2]^2, p) ge 0 then
    u:=1/z^Valuation(L[1]*L[2],p)*z^Valuation(disc,p);
    f:=RR!1;
    if c lt Kmu(p,L,mu) then
      for n in [1..c] do
        if #Lmu(p,L,mu,n) mod 2 eq 0 then
          d:=Integers()!(2*dmu(p,L,mu,n));
          f:=f+(1-1/z^2)*vmu(p,L,mu,n)*z^d*X^n;
        end if;
      end for;
      d:=Integers()!(2*dmu(p,L,mu,c+1));
      f:=f+vmu(p,L,mu,c+1)*fmu(p,L,mu,tmu(p,L,mu,m))*z^d*X^(c+1);
    else
      for n in [1..Kmu(p,L,mu)] do
        if #Lmu(p,L,mu,n) mod 2 eq 0 then
          d:=Integers()!(2*dmu(p,L,mu,n));
          f:=f+(1-1/z^2)*vmu(p,L,mu,n)*z^d*X^n;
        end if;
      end for;
    end if;
    S:=u*f;
  end if;
  for n in [0..Degree(S)] do
    num:=P!Numerator(Coefficient(S,n)) mod P!(z^2-p);
    den:=P!Denominator(Coefficient(S,n)) mod P!(z^2-p);
    T:=T+num/den*X^n;
  end for;
  if disc eq -3 and p eq 3 then
    T:=3*T;
  end if;
  return T;
end function;

dmu2:=function(L,mu,m)
  S:=m;
  for n in Hmu(2,mu) do
    ell:=Valuation(L[n],2);
    S:=S+Min(0,(ell-m+1)/2);
  end for;
  return S;
end function;

epsmu:=function(L,mu,m)
  S:=1;
  for n in Lmu(2,L,mu,m-1) do
    u:=L[n]/2^Valuation(L[n],2);
    S:=S*u;
  end for;
  return S;
end function;

deltamu:=function(L,mu,k)
  S:=1;
  for n in Hmu(2,mu) do
    if Valuation(L[n],2) eq k-1 then
      S:=0;
    end if;
  end for;
  return S;
end function;

Qpmu:=function(L,mu)
  S:=0;
  for n in [i: i in [1..2]| not (i in Hmu(2,mu))] do
    S:=S+L[n]*mu[n]^2;
  end for;
  return S;
end function;

numu:=function(L,mu,m,k)
  S:=(m-Qpmu(L,mu))*2^(3-k);
  for n in Hmu(2,mu) do
    ell:=Valuation(L[n],2);
    eps:=L[n]/2^ell;
    if ell lt k-1 then
      S:=S-eps;
    end if;
  end for;
  return S;
end function;

Kmu2:=function(L,mu)
  S:=1000000000000000;
  for n in [1..2] do
    e:=Valuation(mu[n],2);
    if e eq -1 then
      S:=Min(S,Valuation(L[n],2)+1);
    end if;
    if e lt -1 then
      S:=Min(S,Valuation(L[n],2)+e+1);
    end if;
  end for;
  return S;
end function;

Wmu2:=function(disc,L,mu,m)
  u:=z^(-Valuation(L[1]*L[2],2)-2)*z^Valuation(disc,2);
  if Valuation(m-L[1]*mu[1]^2-L[2]*mu[2]^2, 2) ge 0 then
    S:=1;
  else
    S:=0;
  end if;
  KMax:=Valuation(m-Qpmu(L,mu),2)+3;
  for k in [1..Min(Kmu2(L,mu),KMax)] do
    v:=deltamu(L,mu,k)*z^(Integers()!(2*dmu2(L,mu,k)));
    if #Lmu(2,L,mu,k-1) mod 2 eq 1 then
      v:=v/z^3;
      if Valuation(numu(L,mu,m,k),2) eq 0 then
         v:=v*KroneckerSymbol(8,
           Numerator(numu(L,mu,m,k)*epsmu(L,mu,k))*
           Denominator(numu(L,mu,m,k)*epsmu(L,mu,k)));
      else
        v:=0;
      end if;
    else
      v:=v/z^2*KroneckerSymbol(8,Numerator(epsmu(L,mu,k))*
        Denominator(epsmu(L,mu,k)));
      if Valuation(numu(L,mu,m,k),2) lt 2 then
        v:=0;
      end if;
      if Valuation(numu(L,mu,m,k),2) eq 2 then
        v:=-v;
      end if;
    end if;
    S:=S+v*X^k;
  end for;
  S:=RR!(u*S); T:=RR!0;
  for n in [0..Min(Kmu2(L,mu),KMax)] do
    num:=P!Numerator(Coefficient(S,n)) mod P!(z^2-2);
    den:=P!Denominator(Coefficient(S,n)) mod P!(z^2-2);
    T:=T+num/den*X^n;
  end for;
  return T;
end function;

Kmu22:=function(Sgcd,SS,mu)
  if #Hmu(2,mu) eq 2 then
    return 10000000000000000;
  else
    return Valuation(Sgcd,2)+Min(Valuation(mu[1],2),Valuation(mu[2],2));
  end if;
end function;

dmu22:=function(Sgcd,SS,mu,k)
  if #Hmu(2,mu) eq 2 then
    return Min(Valuation(Sgcd,2),k);
  else
    return k;
  end if;
end function;

pmu:=function(Sgcd,SS,mu,k)
  if #Hmu(2,mu) eq 2 then
    return (-1)^Min(Valuation(Sgcd,2)-k,0);
  else
    return 1;
  end if;
end function;

Qmu:=function(Sgcd,SS,mu)
  return Sgcd*(SS[1,1]*mu[1]^2/2+SS[1,2]*mu[1]*mu[2]+SS[2,2]*mu[2]^2/2);
end function;

Qpmu2:=function(Sgcd,SS,mu)
  if #Hmu(2,mu) eq 2 then
    return 0;
  else
    return Qmu(Sgcd,SS,mu);
  end if;
end function;

numu2:=function(Sgcd,SS,mu,m,k)
  return (m-Qpmu2(Sgcd,SS,mu))*2^(3-k);
end function;

Wmu22:=function(Sgcd,SS,disc,mu,m)
  u:=1/z^(2*Valuation(Sgcd,2))*z^Valuation(disc,2);
  if Valuation(m-Qmu(Sgcd,SS,mu),2) ge 0 then
    S:=1;
  else
    S:=0;
  end if;
  KMax:=Valuation(m-Qpmu2(Sgcd,SS,mu),2)+3;
  for k in [1..Min(KMax,Kmu22(Sgcd,SS,mu))] do
    d:=Integers()!(2*dmu22(Sgcd,SS,mu,k)-2);
    if Valuation(numu2(Sgcd,SS,mu,m,k),2) eq 2 then
      S:=S-pmu(Sgcd,SS,mu,k)*z^d*X^k;
    end if;
    if Valuation(numu2(Sgcd,SS,mu,m,k),2) gt 2 then
      S:=S+pmu(Sgcd,SS,mu,k)*z^d*X^k;
    end if;
  end for;
  S:=RR!(u*S); T:=RR!0;
  for n in [0..Min(Kmu22(Sgcd,SS,mu),KMax)] do
    num:=P!Numerator(Coefficient(S,n)) mod P!(z^2-2);
    den:=P!Denominator(Coefficient(S,n)) mod P!(z^2-2);
    T:=T+num/den*X^n;
  end for;
  return T;
end function;

Wmu23:=function(Sgcd,SS,disc,mu,m)
  u:=1/z^(2*Valuation(Sgcd,2))*z^Valuation(disc,2);
  if Valuation(m-Qmu(Sgcd,SS,mu),2) ge 0 then
    S:=1;
  else
    S:=0;
  end if;
  KMax:=Valuation(m-Qpmu2(Sgcd,SS,mu),2)+3;
  for k in [1..Min(KMax,Kmu22(Sgcd,SS,mu))] do
    d:=Integers()!(2*dmu22(Sgcd,SS,mu,k)-2);
    if Valuation(numu2(Sgcd,SS,mu,m,k),2) eq 2 then
      S:=S-z^d*X^k;
    end if;
    if Valuation(numu2(Sgcd,SS,mu,m,k),2) gt 2 then
      S:=S+z^d*X^k;
    end if;
  end for;
  S:=RR!(u*S); T:=RR!0;
  for n in [0..Min(Kmu22(Sgcd,SS,mu),KMax)] do
    num:=P!Numerator(Coefficient(S,n)) mod P!(z^2-2);
    den:=P!Denominator(Coefficient(S,n)) mod P!(z^2-2);
    T:=T+num/den*X^n;
  end for;
  return T;
end function;

// This is Kappa_0(m)
GetHeight:=function(disc,Lplus,m)
  printf "Kappa_0 of %o\n", m;
  Mplus:=LCoefficients(Lplus);
  M:=Matrix(Rationals(),3,3,[LCoefficients(Proj(i,{Lplus})),
     LCoefficients(Proj(j,{Lplus})),LCoefficients(Proj(k,{Lplus}))]);
  B:=Basis(RowSpace(M));
  LL:=LatticeWithBasis(3,[1,0,0,0,1,0,0,0,1]);
  MM:=LatticeWithBasis(3,[B[1][1],B[1][2],B[1][3],B[2][1],
       B[2][2],B[2][3]]);
  B:=Basis(LL meet MM);
  V1:=B[1]; V2:=B[2];
  w1:=V1[1]*L[1]+V1[2]*L[2]+V1[3]*L[3];
  w2:=V2[1]*L[1]+V2[2]*L[2]+V2[3]*L[3];
  M0:=Matrix(Rationals(),2,2,[QForm(w1,w1),QForm(w1,w2),
     QForm(w2,w1),QForm(w2,w2)]);
  SS,U:=LLLGram(M0);
  SS:=MatrixAlgebra(Integers(),2)!SS;
  Sgcd:=GCD([SS[1,1],SS[1,2],SS[2,1],SS[2,2]]);
  printf "\n\t Sgcd = %o\n", Sgcd;
  SS:=Matrix(Integers(),2,2,[SS[1,1] div Sgcd, SS[1,2] div Sgcd,
      SS[2,1] div Sgcd, SS[2,2] div Sgcd]);
  Mminus1:=U[1,1]*V1+U[1,2]*V2;
  Mminus2:=U[2,1]*V1+U[2,2]*V2;
  Mm:=Matrix(Rationals(),2,3,
    [Mminus1[1],Mminus1[2],Mminus1[3],Mminus2[1],Mminus2[2],Mminus2[3]]);
  Lminus1:=Mminus1[1]*L[1]+Mminus1[2]*L[2]+Mminus1[3]*L[3];
  Lminus2:=Mminus2[1]*L[1]+Mminus2[2]*L[2]+Mminus2[3]*L[3];
  W:=LCoefficients(Lplus);
  M:=Matrix(Integers(),3,3,
     [[Integers()!W[1],Integers()!W[2],Integers()!W[3]],
      [Mminus1[1],Mminus1[2],Mminus1[3]],
      [Mminus2[1],Mminus2[2],Mminus2[3]]]);
  M0,V0:=HermiteForm(M);

  S:=1; pps:={};
  exps:=AssociativeArray(); exps[0]:=0;
  for r in [0..M0[1,1]-1] do
    for s in [0..M0[2,2]-1] do
      for t in [0..M0[3,3]-1] do
        eta:=r*L[1]+s*L[2]+t*L[3];
        etaminus:=Proj(eta,{Lplus});
        vv:=LCoefficients(etaminus);
        nu:=Solution(Mm,Vector(Rationals(),3,[vv[1],vv[2],vv[3]]));
        etaplus:=eta-etaminus;
        u:=Rationals()!(etaplus*Conjugate(Lplus)/Norm(Lplus));
        v:=Sqrt(m/Norm(Lplus));
        for n in [Ceiling(-v-u)..Floor(v-u)] do
          x:=(u+n)*Lplus;
          mp:=m-Norm(x); pp:=0; B:=1;
          printf "\n\t mu_minus = %o, m - Q(x) = %o\n", etaminus, mp;
          if mp eq 0 then
            if IsIntegral(nu[1]) and IsIntegral(nu[2]) then
              if 0 notin pps then
                Include(~pps,0);
                exps[0]:=0;
              end if;
              exps[0]:=exps[0]+1;
            end if;
            dd:=FundamentalDiscriminant(disc);
            if #PrimeDivisors(disc div dd) eq 1 then
              LL:=Factorization(disc div dd);
              p:=LL[1][1]; e:=LL[1][2] div 2;
              if p notin pps then
                Include(~pps,p); exps[p]:=0;
              end if;
              exps[p]:=exps[p]+2*p^(1-e)/(p-KroneckerSymbol(dd,p));
            end if;
            cc:=disc;
            if cc mod 4 eq 0 then cc:=cc div 4; end if;
            x1,y1:=SquarefreeFactorization(cc);
            x2,y2:=SquarefreeFactorization(m);
            if y1 eq 1 then
              DD:=Factorization(y2);
              for L in DD do
                p:=L[1];
                if p notin pps then
                  Include(~pps,p);
                  exps[p]:=0;
                end if;
                exps[p]:=exps[p]+2*L[2];
              end for;
            end if;
          else
            cc:=FundamentalDiscriminant(disc);
            ee,dd:=SquarefreeFactorization(disc div cc);
            
            for p in PrimeDivisors(Sgcd*disc*Numerator(mp)) do
              if not (p in pps) then
                Include(~pps,p);
                exps[p]:=0;
              end if;
              if p gt 2 then
                M,A:=GetDiag(Sgcd,SS,p); LL:=[M[1,1]/2, M[2,2]/2];
                mu:=nu*A^-1;
                f:=Wmu(disc,p,LL,mu,mp);
                if (f mod (X-1) eq 0) and (pp eq 0) then
                  pp:=p;
                  if f eq 0 then d:=0; else d:=Degree(f); end if;
                  B:=B*&+([-Coefficient(f,i)*i: i in [0..d]])/
                     (1-LegendreSymbol(cc,p)/p);
                else
                  B:=B*(f mod (X-1))/(1-LegendreSymbol(cc,p)/p);
                end if;
              else
                M,A:=GetDiag(Sgcd,SS,2);
                if M[1,2] eq 0 then
                  LL:=[M[1,1]/2,M[2,2]/2];
                  mu:=nu*A^-1;
                  f:=Wmu2(disc,LL,mu,mp);
                else
                  if (Valuation(SS[1,1],2) eq 1) and
                     (Valuation(SS[2,2],2) eq 1) then
                    f:=Wmu22(Sgcd,SS,disc,nu,mp);
                  else
                    f:=Wmu23(Sgcd,SS,disc,nu,mp);
                  end if;
                end if;
                if (f mod (X-1) eq 0) and (pp eq 0) then
                  pp:=2;
                  if f eq 0 then d:=0; else d:=Degree(f); end if;
                  B:=B*&+([-Coefficient(f,i)*i: i in [0..d]])/
                     (1-KroneckerSymbol(cc,2)/2);
                else
                  B:=B*(f mod (X-1))/(1-KroneckerSymbol(cc,2)/2);
                end if;
              end if;
            end for;
            cc:=FundamentalDiscriminant(disc);
            ee,dd:=SquarefreeFactorization(disc div cc);
            if cc eq -4 then B:=B*2; end if;
            if cc eq -3 then B:=B*3; end if;
            printf "\t adding %m log %m\n", Rationals()!(B*2/ClassNumber(cc)/dd), pp;
            exps[pp]:=exps[pp]+B*2/ClassNumber(cc)/dd; // This is adding the part from Kappa_eta(mu)
          end if;
        end for;
      end for;
    end for;
  end for;
  return pps, exps;
end function;

GetHeight1:=function(disc,Lplus,m)
  printf "Kappa_eta of %o\n", m;
  Mplus:=LCoefficients(Lplus);
  M:=Matrix(Rationals(),3,3,[LCoefficients(Proj(i,{Lplus})),
     LCoefficients(Proj(j,{Lplus})),LCoefficients(Proj(k,{Lplus}))]);
  B:=Basis(RowSpace(M));
  LL:=LatticeWithBasis(3,[1,0,0,0,1,0,0,0,1]);
  MM:=LatticeWithBasis(3,[B[1][1],B[1][2],B[1][3],B[2][1],
       B[2][2],B[2][3]]);
  B:=Basis(LL meet MM);
  V1:=B[1]; V2:=B[2];
  w1:=V1[1]*L[1]+V1[2]*L[2]+V1[3]*L[3];
  w2:=V2[1]*L[1]+V2[2]*L[2]+V2[3]*L[3];
  M0:=Matrix(Rationals(),2,2,[QForm(w1,w1),QForm(w1,w2),
     QForm(w2,w1),QForm(w2,w2)]);
  SS,U:=LLLGram(M0);
  SS:=MatrixAlgebra(Integers(),2)!SS;
  Sgcd:=GCD([SS[1,1],SS[1,2],SS[2,1],SS[2,2]]);
  SS:=Matrix(Integers(),2,2,[SS[1,1] div Sgcd, SS[1,2] div Sgcd,
      SS[2,1] div Sgcd, SS[2,2] div Sgcd]);
  Mminus1:=U[1,1]*V1+U[1,2]*V2;
  Mminus2:=U[2,1]*V1+U[2,2]*V2;
  Mm:=Matrix(Rationals(),2,3,
    [Mminus1[1],Mminus1[2],Mminus1[3],Mminus2[1],Mminus2[2],Mminus2[3]]);
  Lminus1:=Mminus1[1]*L[1]+Mminus1[2]*L[2]+Mminus1[3]*L[3];
  Lminus2:=Mminus2[1]*L[1]+Mminus2[2]*L[2]+Mminus2[3]*L[3];
  W:=LCoefficients(Lplus);
  M:=Matrix(Integers(),3,3,
     [[Integers()!W[1],Integers()!W[2],Integers()!W[3]],
      [Mminus1[1],Mminus1[2],Mminus1[3]],
      [Mminus2[1],Mminus2[2],Mminus2[3]]]);
  M0,V0:=HermiteForm(M);
  eta1:=L[1]; eta2:=L[2]; eta3:=L[3];

  S:=1; pps:={}; tt:=false;
  exps:=AssociativeArray(); exps[0]:=0;
if Integers()!(4*m) mod 4 eq 0 or Integers()!(4*m) mod 4 eq 3 then
  for r in [0..M0[1,1]-1] do
    for s in [0..M0[2,2]-1] do
      for t in [0..M0[3,3]-1] do
        if Integers()!(4*m) mod 4 eq 3 then
          eta:=r*eta1+s*eta2+t*eta3+D*DualL[3];
        else
          eta:=r*eta1+s*eta2+t*eta3;
        end if;
        etaminus:=Proj(eta,{Lplus});
        vv:=LCoefficients(etaminus);
        nu:=Solution(Mm,Vector(Rationals(),3,[vv[1],vv[2],vv[3]]));
        etaplus:=eta-etaminus;
        u:=Rationals()!(etaplus*Conjugate(Lplus)/Norm(Lplus));
        v:=Sqrt(m/Norm(Lplus));
        for n in [Ceiling(-v-u)..Floor(v-u)] do
          x:=(u+n)*Lplus;
          mp:=m-Norm(x); pp:=0; B:=1;
          printf "\tm - Q(x) = %o\n", mp;
          if mp eq 0 then
            if IsIntegral(nu[1]) and IsIntegral(nu[2]) then
              if 0 notin pps then
                Include(~pps,0);
                exps[0]:=0;
              end if;
              exps[0]:=exps[0]+1;
            end if;
            tt:=true;
            dd:=FundamentalDiscriminant(disc);
            if #PrimeDivisors(disc div dd) eq 1 then
              LL:=Factorization(disc div dd);
              p:=LL[1][1]; e:=LL[1][2] div 2;
              if p notin pps then
                Include(~pps,p); exps[p]:=0;
              end if;
              exps[p]:=exps[p]+0*p^(1-e)/(p-KroneckerSymbol(dd,p));
            end if;
          else
            for p in PrimeDivisors(2*Sgcd*disc*Numerator(mp)) do
              if not (p in pps) then
                Include(~pps,p);
                exps[p]:=0;
              end if;
              if p gt 2 then
                M,A:=GetDiag(Sgcd,SS,p); LL:=[M[1,1]/2, M[2,2]/2];
                mu:=nu*A^-1;
                f:=Wmu(disc,p,LL,mu,mp);
                if (f mod (X-1) eq 0) and (pp eq 0) then
                  pp:=p;
                  if f eq 0 then d:=0; else d:=Degree(f); end if;
                  B:=B*&+([-Coefficient(f,i)*i: i in [0..d]])/
                     (1-LegendreSymbol(disc,p)/p);
                else
                  B:=B*(f mod (X-1))/(1-LegendreSymbol(disc,p)/p);
                end if;
              else
                M,A:=GetDiag(Sgcd,SS,2);
                if M[1,2] eq 0 then
                  LL:=[M[1,1]/2,M[2,2]/2];
                  mu:=nu*A^-1;
                  f:=Wmu2(disc,LL,mu,mp);
                else
                  if (Valuation(SS[1,1],2) eq 1) and
                     (Valuation(SS[2,2],2) eq 1) then
                    f:=Wmu22(Sgcd,SS,disc,nu,mp);
                  else
                    f:=Wmu23(Sgcd,SS,disc,nu,mp);
                  end if;
                end if;
                if (f mod (X-1) eq 0) and (pp eq 0) then
                  pp:=2;
                  if f eq 0 then d:=0; else d:=Degree(f); end if;
                  B:=B*&+([-Coefficient(f,i)*i: i in [0..d]])/
                     (1-KroneckerSymbol(disc,2)/2);
                else
                  B:=B*(f mod (X-1))/(1-KroneckerSymbol(disc,2)/2);
                end if;
              end if;
            end for;
            printf "\tadding %o Log %o\n", B*2/ClassNumber(disc), pp;
            exps[pp]:=exps[pp]+B*2/ClassNumber(disc);
          end if;
        end for;
      end for;
    end for;
  end for;
  if tt then
    dd:=FundamentalDiscriminant(disc);
    if #PrimeDivisors(disc div dd) eq 1 then
      LL:=Factorization(disc div dd);
      p:=LL[1][1]; e:=LL[1][2] div 2;
      if p notin pps then
        Include(~pps,p); exps[p]:=0;
      end if;
      printf "\tadding %o Log %o\n",4*p^(1-e)/(p-KroneckerSymbol(dd,p)), p;
      exps[p]:=exps[p]+4*p^(1-e)/(p-KroneckerSymbol(dd,p));
    end if;
  end if;
end if;
  return pps, exps;
end function;

for r in [0..NMax] do
  for s in [-NMax..NMax] do
    for t in [c: c in [-NMax..NMax]|GCD([r,s,c]) eq 1] do
      beta:=r*L[1]+s*L[2]+t*L[3];
      if Norm(beta) gt 0 then
        d:=Integers()!(-Norm(beta)); u:=4*d;
        if d mod 4 eq 1 then
          if (1+beta)/2 in OO then u:=d; end if;
        end if;
        if (u notin Discs) and u lt 0 then
          Include(~Discs,u);
          Optimals[u]:=beta;
          ClassNos[u]:=ClassNumber(u)*
            &*([1-KroneckerSymbol(u,p): p in PrimeDivisors(D)]);
          if D mod d eq 0 then
            ClassNos[u]:=Integers()!(ClassNos[u]/2);
          else
            ClassNos[u]:=Integers()!(ClassNos[u]/4);
          end if;
          if u eq -4 then ClassNos[-4]:=1; end if;
        end if;
      end if;
    end for;
  end for;
end for;

SingularModuli:=procedure(disc)
  print "====================================================";
  print "The discriminant is",disc,"with class number",ClassNos[disc],
    "and embedding",Optimals[disc];

  S1:=&join([Seqset(L): L in divs1]);
  exps1:=AssociativeArray();
  pps1:=AssociativeArray();
  for n in S1 do
    exps1[n]:=AssociativeArray();
    pps1[n]:={};
  end for;
  for m in S1 do
    pp,ex:=GetHeight1(disc,Optimals[disc],m);
    for p in pp do
      if p notin pps1[m] then
        exps1[m][p]:=ex[p];
        Include(~pps1[m],p);
      else
        exps1[m][p]:=exps1[m][p]+ex[p];
      end if;
    end for;
  end for;

  S:=&join([Seqset(L): L in divs]);
  exps:=AssociativeArray();
  pps:=AssociativeArray();
  for n in S do
    exps[n]:=AssociativeArray();
    pps[n]:={};
  end for;
  for m in S do
    pp,ex:=GetHeight(disc,Optimals[disc],m);
    for p in pp do
      if p notin pps[m] then
        exps[m][p]:=ex[p];
        Include(~pps[m],p);
      else
        exps[m][p]:=exps[m][p]+ex[p];
      end if;
    end for;
  end for;

  for n in [1..#divs] do
    ps:={p: p in offs[n]};
    exs:=AssociativeArray();
    for i in [1..#offs[n]] do
      exs[offs[n][i]]:=offsets[n][i]*4;
    end for;

    for i in [1..#divs[n]] do
      m:=divs[n][i]; wt:=wts[n][i];
      for p in pps[m] do
        if p notin ps then
          exs[p]:=wt*exps[m][p];
          Include(~ps,p);
        else
          exs[p]:=exs[p]+wt*exps[m][p];
        end if;
      end for;
    end for;

    for i in [1..#divs1[n]] do
      m:=divs1[n][i]; wt:=wts1[n][i];
      for p in pps1[m] do
        if p notin ps then
          exs[p]:=wt*exps1[m][p];
          Include(~ps,p);
        else
          exs[p]:=exs[p]+wt*exps1[m][p];
        end if;
      end for;
    end for;

    S:=1;
    if 0 in ps then
      if exs[0] gt 0 then
        S:=0;
        print "For the divisor",divs[n],"with weight",wts[n], S;
      end if;
      if exs[0] lt 0 then
        S:=Infinity();
        print "For the divisor",divs[n],"with weight",wts[n], S;
      end if;
      if exs[0] eq 0 then
        L1:=[]; L2:=[];
        for p in ps do
          if Rationals()!exs[p] gt 0 then
            Append(~L1,<p,ClassNos[disc]*exs[p]/4>);
          end if;
          if Rationals()!exs[p] lt 0 then
            Append(~L2,<p,-ClassNos[disc]*exs[p]/4>);
          end if;
        end for;
        print "For the divisor",divs[n],"with weight",wts[n],
          "the numerator is";
        print L1;
        print "the denominator is";
        print L2;
      end if;
    else
      L1:=[]; L2:=[];
      for p in ps do
        if Rationals()!exs[p] gt 0 then
          Append(~L1,<p,ClassNos[disc]*exs[p]/4>);
        end if;
        if Rationals()!exs[p] lt 0 then
          Append(~L2,<p,-ClassNos[disc]*exs[p]/4>);
        end if;
      end for;
      print "For the divisor",divs[n],"with weight",wts[n],
        "the numerator is";
      print L1;
      print "the denominator is";
      print L2;
    end if;
    if n ne #divs then
      print "--------------------------------------------------";
    end if;
  end for;
end procedure;

for disc in [r: r in [1..2000]| (-r in Discs and ClassNos[-r] le 1)] do
  SingularModuli(-disc);
end for;

