# Fibre-product entries for the level-924 bases (banked 2026-10-11; now PR #82)

**PR #82** (`level-924-fibre-products`, based on #72's branch with #80 merged in) carries the 56
level-924 entries appended to their model files. `recheck_crv.m` is the independent check run for
it: every entry re-read from its stored strings, compared with the trace formula at the 4th and 5th
good primes (112/112). The 22_7 and 26_5 appends are still only here.

`tools/writefp.m` from the branch `fibre-product-models` (PR #72, at commit 4bd3d86), run over the
model files of PR #80 (`models_462_1.m`, `models_6_77.m`, `models_22_21.m` copied into that branch's
`data/models/`), wrote these `.append` files: one `<genus, "CRV", [equations]>` entry per cover that
had no equation but whose index-2 Atkin-Lehner double covers do. The writer keeps an entry only if
the genus of the compositum equals the Shimura-curve genus and its point counts over F_p and F_{p^2}
agree with the Eichler-Selberg trace formula at up to three good primes; every rejection is in
`writefp.log`.

| file | entries | what they fill |
|---|---|---|
| `models_22_21.append` | 28 | **every** empty key of X_0^22(21) (genera 3 to 7) |
| `models_6_77.append` | 28 | **every** empty key of X_0^6(77) (genera 3 to 9) |
| `models_22_7.append` | 1 | one genus-3 cover of X_0^22(7), new since #72 (its double covers came with #70) |
| `models_26_5.append` | 1 | one genus-4 cover of X_0^26(5), likewise |

X_0^462(1) got no entry: its 11 empty keys have no committed double covers to multiply.

Independent check done here: for both level-924 bases, the appended keys are exactly the empty keys
of the committed file, and every stored genus equals the genus in the curve table (0 mismatches).

**Where they go.** The `"CRV"` entry shape and the stage that builds and tests it live on
`fibre-product-stage` (PR #71), not on `main`, so these cannot join PR #80 (which targets `main`).
Once #71 and #80 are merged, append each file to its `data/models/models_D_N.m` and add the
provenance lines; the 22_7 and 26_5 entries belong with PR #72's data. The appends are plain Magma
assignments and can be concatenated as they are.
