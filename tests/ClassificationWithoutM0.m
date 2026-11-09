// The classification uses no CM value at level N > 1, where Schofer's formula has an m = 0 term, and no
// model built from such values.  The two places where such data could enter a verdict are the CM tables
// of the special-fibre test (special_fiber_cm.m) and the list of proved models that the canonical-quadrics
// stage reads (data/models/proved_models.m); the test fails if either holds an entry at N > 1.  A run of
// the whole pipeline with Schofer's formula at N > 1 and the m = 0 multipliers made to raise an error
// reproduced every verdict (docs/CLASSIFICATION_PIPELINE.md, "Independence from the m = 0 term").

import "special_fiber_cm.m" : CM_BASE_DATA;

assert forall(k){k : k in Keys(CM_BASE_DATA) | k[2] eq 1};
proved := eval (Read("data/models/proved_models.m") cat "\nreturn proved_models;");
assert forall(e){e : e in proved | e[2] eq 1};
