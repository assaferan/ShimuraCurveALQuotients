# Probes behind PR #84: the class of complex conjugation on CM points, read off the Eichler order

Sachi's example (2026-10-11): on X_0(15,1)/<w_15> the CM points by the order of discriminant -48
are defined over Q(i); `FieldsOfDefinitionOfCMPointFast` returned Q(i) and Q(sqrt 3).

* `lambda.m` -- `magma -b D:=15 N:=1 d:=-48 m:=5 lambda.m`: for every optimal embedding of the
  order of discriminant d into the Eichler order found in a coefficient box, the norm form of the
  lattice of anticommuting elements divided by -m, its discriminant and reduced class.  At
  (15, 1, -48) all 12 embeddings give the principal class of discriminant -48; at (6, 1, -75) all 24
  give one class (content 2 there, so read the class, not the content).  The box misses embeddings
  with large coefficients, so "0 embeddings" means nothing; the committed intrinsic uses FindLambda.
* `gydriver.m`, `gydriver.log` -- the scratch package's field (same construction as the committed
  `CMConjugationIdeal`) against every Guo-Yang field of degree >= 2 in `tests/CMPoints.m`: 460
  points, the quaternion-algebra condition ambiguous at 12, the lattice class matching Guo-Yang at
  all 12 and the old unique answer at the other 448.  The one recorded disagreement with Guo-Yang,
  -132 on X_0^206(1)/<w_206>, is not an ambiguous case and is unchanged.

Two traps met on the way: the order of discriminant d is Z[sqrt(d/4)] when 4 | d, so searching for
x^2 = d embeds the order of discriminant 4d; and the pipeline's own CM value at -48 carries the
2^(4/3) artefact of the non-fundamental m = 0 term, so the model files cannot arbitrate this case.
