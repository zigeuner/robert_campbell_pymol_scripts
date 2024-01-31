#! /usr/bin/env python
# Copyright (c) 2004 Robert L. Campbell

import sys
from cctbx import sgtbx

def get_symm(sg):

  try:
    space_group_info = sgtbx.space_group_info(symbol=sg)
  except RuntimeError, err:
    print err
    sys.exit(1)

  space_group_info.show_summary()
  print len(space_group_info.group()), space_group_info.group().n_smx()
#  print space_group_info.Info().BuildLookupSymbol()

  # print space_group symmetry information as x,y,z and as matrix row by row
  i = 1
  for smx in space_group_info.group():
    tlist = []
    rlist = []
    r, t = smx.r(), smx.t()
    for j in range(len(r.num())):
      rlist.append(float(r.num()[j]))
    for j in range(len(t.num())):
      tlist.append(t.num()[j]/12.)

    print "%3d %-24s %-50s %-20s" % (i, str(smx), str(rlist), str(tlist))
    i += 1


if (__name__ == "__main__"):
  import sys
  # if no argument given, then list the whole shebang!
  if (len(sys.argv) == 1):
    for i in xrange(230):
      get_symm(i + 1)
  # if an argument is given, then print the information for just that space group
  else:
    for symbol in sys.argv[1:]:
      get_symm(symbol)

