#! /usr/bin/env python
# Copyright (c) 2004 Ralf W. Gross-Kunstleve, Robert L. Campbell
# List all axes in the unit cell.

# usage:
#   python all_axes.py     - show axes for the 230 reference settings.
#   python all_axes.py P2  - show axes for (e.g.) space group P2

# RWGK = Ralf W. Grosse-Kunstleve
# RWGK Some further refinement is required:
# RWGK   - List only the axes of highest order (e.g. only 4, not 4 and 2).
# RWGK   - List only the axes with the smallest intrinsic component
# RWGK     (e.g. list only 3(1), not both 3(1) and 3(2)).
# RWGK See also: comment regarding shift_range below.
# modified a little by rlc

from cctbx import sgtbx
from math import sqrt

import string, re

def str_ev(EV):
  return "[%d,%d,%d]" % EV

def RTMxAnalysis(M):
  r_info = sgtbx.rot_mx_info(M.r())
  t_info = sgtbx.translation_part_info(M)
  t_intrinsic = str(t_info.intrinsic_part().mod_positive())
  t_intrinsic_float = t_info.intrinsic_part().mod_positive().as_double()
  t_shift = str(t_info.origin_shift().mod_positive())
#  t_location = str(t_info.location_part().mod_positive())
#dbg#  print "t_intrinsic,t_shift",t_intrinsic,t_shift
  if (r_info.type() == 1):
    return ("1", "-", "-", "-")
  elif (r_info.type() == -1):
    return (str(r_info.type()), "-", "-", "(%s)" % (t_shift,))
  elif (abs(r_info.type()) == 2):
    return (str(r_info.type()), str_ev(r_info.ev()), "(%s)" % (t_intrinsic,), "(%s)" % (t_shift,))
#    return (str(r_info.type()), str_ev(r_info.ev()), "(%s)" % (t_intrinsic,), "(%s)" % (t_shift,), "[%s]" % (t_location,))
#
# restrict to positive sense, if desired:
#  elif (r_info.sense() > 0):
# or not ...
  else:
    trans = sqrt(reduce(lambda x,y:x+y,(map(lambda x,y:(y-x)*(y-x),(0,0,0),t_intrinsic_float))))
    if trans == 0:
      return (str(r_info.type()), str_ev(r_info.ev()), "(%s)" % (t_intrinsic,), "(%s)" % (t_shift,))
      #return (str(r_info.type()), str_ev(r_info.ev()), "(%s)" % (t_intrinsic,), "(%s)" % (t_shift,), "[%s]" % (t_location,))
    else:
      return (str(r_info.type())+"^"+str(r_info.sense()), str_ev(r_info.ev()), "(%s)" % (t_intrinsic,), "(%s)" % (t_shift,))
      #return (str(r_info.type())+"^"+str(r_info.sense()), str_ev(r_info.ev()), "(%s)" % (t_intrinsic,), "(%s)" % (t_shift,), "[%s]" % (t_location,))

def list_all_axes(space_group_symbol=None, space_group_info=None):
  assert space_group_symbol is None or space_group_info is None
  shift_range = 1 # RWGK Works for the 230 reference settings; it is not
                  # RWGK clear to me (rwgk) what value is needed in general.
  if (space_group_symbol is not None):
    space_group_info = sgtbx.space_group_info(symbol=space_group_symbol)
  space_group_info.show_summary()
  print ""
#  print 'Hall symbol: ', space_group_info.type().hall_symbol()
  print "Rotation type, Axis direction, Intrinsic part, Origin shift"
#dbg#  print "Rotation type, Axis direction, Intrinsic part, Origin shift, Location part"
  axes_dict = {}
  for smx in space_group_info.group():
    r = smx.r()
    t = smx.t()
    shift = [0,0,0]
    for shift[0] in range(-shift_range,shift_range+1):
      for shift[1] in range(-shift_range,shift_range+1):
        for shift[2] in range(-shift_range,shift_range+1):
          ts = t.plus(sgtbx.tr_vec(shift, 1)).new_denominator(t.den())
          m = sgtbx.rt_mx(r, ts)
          rtmxanal = RTMxAnalysis(m)
          axes_dict[rtmxanal] = 0
# print operator and rtmx analysis to check all operators if desired
#dbg#          print m, rtmxanal
  axes_list = axes_dict.keys()
  axes_list.sort()
  for a in axes_list:
    print a
  #return axes_dict


if (__name__ == "__main__"):
  import sys
  if (len(sys.argv) == 1):
    for i in range(230):
      list_all_axes(i + 1)
  else:
    for symbol in sys.argv[1:]:
      list_all_axes(symbol)
