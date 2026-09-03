#!/usr/bin/env python3
"""Parse the 42 'CM-values of X_0^D(N)' appendix tables out of the Guo-Yang arXiv
source (1510.06193) into structured JSON. Primary-column (s) values are always plain
rationals or infinity in every table sampled; secondary columns are kept as raw LaTeX
strings for now (mostly radicals), with a best-effort rational parse where possible.
"""
import re
import json
import sys

SRC = "ShimuraCurves-arxiv.tex"

with open(SRC) as f:
    text = f.read()

# Tables are NOT one-per-\begin{table} float (several small ones share a float), and the
# \bigskip/\caption ordering after \end{tabular} is inconsistent across the document. So:
# find every \begin{tabular}...\end{tabular} span (they don't nest, so non-greedy .*? is
# exact), find every CM-values caption's position, and pair each caption with the nearest
# tabular ending before it.
tabular_spans = [m.span() for m in re.finditer(r"\\begin\{tabular\}.*?\\end\{tabular\}", text, re.S)]
caption_re = re.compile(r"\\caption\{CM-values of \$X\^\{(\d+)\}_0\((\d+)\)\$\}")

table_blocks = []
for cm in caption_re.finditer(text):
    cstart = cm.start()
    best = None
    for (s0, e0) in tabular_spans:
        if e0 <= cstart and (best is None or e0 > best[1]):
            best = (s0, e0)
    assert best is not None, f"no tabular found before caption at {cstart}"
    s0, e0 = best
    table_blocks.append(text[s0:cm.end()])

def clean_math(s):
    s = s.strip()
    s = re.sub(r"\\ ", " ", s)
    s = re.sub(r"\s+", " ", s)
    return s.strip()

RAT_RE = re.compile(
    r"^\$?\{?\s*(-?\d+)\s*(?:/\s*\{?(\d+)\}?)?\s*\}?\$?$"
)

def parse_scalar(cell):
    """Try to parse a table cell as an exact rational or 'oo'. Returns (kind, value)
    kind in {'rat','oo','other'}; value is a python Fraction numerator/denominator
    tuple for 'rat', None otherwise."""
    c = cell.strip()
    c = c.replace("$", "").replace("{", "").replace("}", "").strip()
    c = re.sub(r"\s+", "", c)
    if c in (r"\infty",):
        return ("oo", None)
    m = re.match(r"^(-?\d+)(?:/(-?\d+))?$", c)
    if m:
        num = int(m.group(1))
        den = int(m.group(2)) if m.group(2) else 1
        return ("rat", [num, den])
    return ("other", c)

results = []
for block in table_blocks:
    cap = re.search(r"\\caption\{CM-values of \$X\^\{(\d+)\}_0\((\d+)\)\$\}", block)
    if not cap:
        continue
    D, N = int(cap.group(1)), int(cap.group(2))

    tab = re.search(r"\\begin\{tabular\}.*?\\end\{tabular\}", block, re.S)
    body = tab.group(0)
    body = body.split("\\begin{tabular}")[1]
    body = body.split("\\end{tabular}")[0]
    # drop the {|c|c|...} spec line (first line before first \hline or content)
    body = re.sub(r"^\{[^}]*\}", "", body.strip())

    # join into logical rows: split on '\\' (row terminator), each row's raw text may
    # itself contain embedded newlines/wrapped continuations -- that's fine, we already
    # have the full row text between one '\\' and the next.
    raw_rows = body.split("\\\\")
    rows_text = []
    for r in raw_rows:
        r = r.replace("\\hline", "")
        r = r.strip()
        if not r:
            continue
        rows_text.append(r)

    if len(rows_text) < 3:
        results.append({"D": D, "N": N, "error": "too few rows", "nrows": len(rows_text)})
        continue

    header1 = rows_text[0]
    header2 = rows_text[1]
    ncols_guess = header1.count("&") + 1

    def split_cols(r):
        return [clean_math(x) for x in r.split("&")]

    h1cols = split_cols(header1)
    h2cols = split_cols(header2)

    data_rows = []
    parse_issues = []
    for r in rows_text[2:]:
        cols = split_cols(r)
        if len(cols) != ncols_guess:
            parse_issues.append({"raw": r, "ncols": len(cols)})
            continue
        disc_cell = cols[0]
        dm = re.search(r"-?\d+", disc_cell.replace("{", "").replace("}", ""))
        if not dm:
            parse_issues.append({"raw": r, "reason": "no disc"})
            continue
        disc = int(dm.group(0))
        vals = []
        for cell in cols[1:]:
            kind, val = parse_scalar(cell)
            vals.append({"kind": kind, "val": val, "raw": cell})
        data_rows.append({"disc": disc, "vals": vals})

    results.append({
        "D": D, "N": N,
        "header1": h1cols, "header2": h2cols,
        "ncols": ncols_guess - 1,
        "rows": data_rows,
        "parse_issues": parse_issues,
    })

with open("tables.json", "w") as f:
    json.dump(results, f, indent=1)

# Summary
print(f"{'D':>4} {'N':>3} {'ncols':>5} {'nrows':>5} {'issues':>6}  {'rational-only col0 count'}")
tot_issues = 0
for t in results:
    if "error" in t:
        print(t)
        continue
    n_rat0 = sum(1 for r in t["rows"] if r["vals"][0]["kind"] in ("rat", "oo"))
    print(f"{t['D']:>4} {t['N']:>3} {t['ncols']:>5} {len(t['rows']):>5} {len(t['parse_issues']):>6}  {n_rat0}/{len(t['rows'])}")
    tot_issues += len(t["parse_issues"])
print(f"\n{len(results)} tables parsed, {tot_issues} total row parse issues")
