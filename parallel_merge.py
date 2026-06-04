#!/usr/bin/env python3
# Merges chunk files produced by parallel_filter_worker.m into a single dat file.
#
# Usage: python3 parallel_merge.py <chunks_dir> <total_chunks> <output_dat>
#
# Each chunk file is a Magma sequence literal:
#   [ PowerStructure(ShimuraQuot) | CreateShimuraQuot(...), ..., CreateShimuraQuot(...) ]
#
# We strip the closing ] from all but the last chunk, and strip the
# opening header from all but the first chunk, joining them with commas.

import sys

chunks_dir   = sys.argv[1]
total_chunks = int(sys.argv[2])
output_dat   = sys.argv[3]

chunks = []
for c in range(1, total_chunks + 1):
    path = f"{chunks_dir}/chunk_{c}_of_{total_chunks}.dat"
    chunks.append(open(path).read())
    print(f"  read chunk {c}: {path}")

# First chunk: keep as-is except strip the trailing ']'
first = chunks[0]
last_bracket = first.rfind(']')
body = first[:last_bracket].rstrip('\n')  # ends with *])

# Middle + last chunks: skip the header line, strip trailing ']' (except for the last)
for i, chunk in enumerate(chunks[1:], start=1):
    after_header = chunk.split('\n', 1)[1]   # drop "[ PowerStructure(ShimuraQuot) |"
    if i < total_chunks - 1:
        # Not the final chunk: strip its trailing ']' too
        last_bracket = after_header.rfind(']')
        body += ',\n' + after_header[:last_bracket].rstrip('\n')
    else:
        # Final chunk: keep the trailing ']'
        body += ',\n' + after_header

with open(output_dat, 'w') as f:
    f.write(body)

print(f"Wrote merged result to {output_dat}")
