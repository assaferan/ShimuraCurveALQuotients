#!/usr/bin/env python3
"""Compare the per-support-class multipliers of the production run (M0MultipliersBySupport lines at
verbosity 2 in a genmodels log) with those of exactm0.m (CMP:=0) on the same base, as multisets of
per-form tuples (the production log lists the forms in the order SchoferFormula first met them).

    python3 cmpprod.py <genmodels log> <exactm0 log>
"""
import re, sys
from collections import Counter

prod_re = re.compile(r"\(1/2\) c_eta\(0\) at the \d+ cosets supported at \{([^}]*)\}: (-?\d+(?:/\d+)?)")
exact_re = re.compile(r"form \d+: \[(.*?)\]\s+\(", re.S)
class_re = re.compile(r"<\{([^}]*)\}, (-?\d+(?:/\d+)?)>")

def key(s):
    return tuple(sorted(int(p) for p in s.replace(" ", "").split(",") if p))

prod = open(sys.argv[1]).read()
items = prod_re.findall(prod)
classes = sorted({key(c) for c, _ in items if key(c)})
k = len(classes)
assert len(items) % k == 0, (len(items), k)
prod_forms = Counter()
for i in range(0, len(items), k):
    chunk = items[i:i + k]
    d = {key(c): v for c, v in chunk}
    assert len(d) == k, chunk
    prod_forms[tuple(d[c] for c in classes)] += 1

exact = open(sys.argv[2]).read()
exact_forms = Counter()
for body in exact_re.findall(exact):
    d = {key(c): v for c, v in class_re.findall(body) if key(c)}
    exact_forms[tuple(d[c] for c in classes)] += 1

print("classes", classes)
print("production forms", sum(prod_forms.values()), "exact forms", sum(exact_forms.values()))
only_prod = prod_forms - exact_forms
only_exact = exact_forms - prod_forms
print("in production only:", dict(only_prod))
print("in exact only:", dict(only_exact))
print("AGREE" if not only_prod and not only_exact else "DIFFER")
