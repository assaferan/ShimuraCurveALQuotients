# Models of subhyperelliptic Atkin–Lehner quotients

Explicit equations for the sub-hyperelliptic (genus ≤ 1, or hyperelliptic) Atkin–Lehner
quotients `X_0(D,N)/W`, computed via the Guo–Yang Borcherds-form method
(`AllEquationsAboveCovers` in `EquationsCovers.m`), generated with polymake (run through
`sage -sh`) for the form-ring solver.

## Files

One file per star base: `models_D_N.m` holds the models for every sub-hyperelliptic cover of
`X_0(D,N)*`. Each is a self-contained, Magma-`load`-able script:

```magma
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
// key = Sort(W) (the Atkin–Lehner subgroup, as a sorted integer list);
// value = list of <genus, f, h>, the model  y^2 + h*y = f(x)   (h is usually 0).
models[[Integers()|1]]   := [* <3, P![...], P![]> *];
models[[Integers()|1,2]] := [* <2, P![...], P![]> *];
```

A list of several `<genus,f,h>` for one key means several model representations were found
(e.g. multiple conic forms for a genus-0 cover); they are isomorphic. Non-hyperelliptic
covers (pointless conics etc.) are stored as `<genus, "CRV", [defining-poly strings]>`.

To use one:
```magma
load "data/models/models_6_11.m";
for k in Keys(models) do for m in models[k] do
  g, f := Explode(<m[1], m[2]>);   // genus g, hyperelliptic polynomial f
end for; end for;
```

## Scope

Generated for the tractable star bases: **D ≥ 2, squarefree N, form-ring level (M = 4·D0,
D0 = D·N / 2^v₂(D)) with ≤ 16 divisors.** D = 1 (classical modular curves) is not handled by
this quaternionic method; larger levels (24+ divisors) are behind a polymake-speed wall.
Sparse/degenerate bases (too few CM points) fail the cover-equation step and are absent.

Verified against Guo–Yang's published tables (e.g. X_0(6,11)) via `tests/X0_*.m`.
