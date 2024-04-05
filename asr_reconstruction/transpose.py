#!/usr/bin/env python
import sys
import os
import re

try:
  in_fn = sys.argv[1]
except ValueError as ve:
  sys.exit(-1)

sz = os.path.getsize(in_fn)
fps = [ open(in_fn, 'r') ]
bufsz = 1000 * 2**20 # MB

with open(in_fn, 'r') as f:
  offset = 0
  buf = f.read(bufsz)
  while len(buf)>0:
    for p in [m.start() for m in re.finditer(r'\n', buf)]:
      new_pos = offset+p+1
      if new_pos < sz:
        new_fp = open(in_fn, 'r')
        new_fp.seek(new_pos, 0)
        fps.append(new_fp)

    offset = f.tell()
    buf = f.read(bufsz)

while sz > 0:
  for fi, f in enumerate(fps):
    byte = f.read(1)
    sz -= 1
    while byte:
      if byte == '\t' or byte == '\n':
        if fi != len(fps) - 1: sys.stdout.write('\t')
        break
      sys.stdout.write('%s' % (byte))
      byte = f.read(1)
      sz -= 1
  sys.stdout.write('\n')

for f in fps:
  f.close()