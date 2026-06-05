#!/usr/bin/env python3
# Merges chunk files produced by parallel_filter_worker.m into a single dat file.
#
# Usage: python3 parallel_merge.py <chunks_dir> <total_chunks> <output_dat>
#
# Each chunk file is a Magma sequence literal:
#   [ PowerStructure(ShimuraQuot) | CreateShimuraQuot(...), ..., CreateShimuraQuot(...) ]
# or empty:
#   [ PowerStructure(ShimuraQuot) | ]

import sys

chunks_dir   = sys.argv[1]
total_chunks = int(sys.argv[2])
output_dat   = sys.argv[3]

chunks = []
for c in range(1, total_chunks + 1):
    path = f"{chunks_dir}/chunk_{c}_of_{total_chunks}.dat"
    chunks.append(open(path).read())
    print(f"  read chunk {c}: {path}")

def extract_elements(chunk):
    """Return the elements string (between '| ' and the final ']'), or '' if empty."""
    pipe = chunk.find('| ')
    last = chunk.rfind(']')
    if pipe == -1 or last == -1:
        return ''
    return chunk[pipe + 2:last].strip()

header = chunks[0].split('\n', 1)[0]  # "[ PowerStructure(ShimuraQuot) |"

parts = []
for c in chunks:
    e = extract_elements(c)
    if e:
        parts.append(e)

if parts:
    body = header + '\n' + ',\n'.join(parts) + ' ]'
else:
    body = header + ' ]'

with open(output_dat, 'w') as f:
    f.write(body)

print(f"Wrote merged result to {output_dat}")
