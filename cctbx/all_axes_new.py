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

from cctbx import sgtbx
#from cctbx.misc.python_utils import list_plus

import numpy as N
import string, re

def list_plus(lhs, rhs):
  return [l + r for l, r in zip(lhs, rhs)]

def list_minus(lhs, rhs):
  return [l - r for l, r in zip(lhs, rhs)]

def list_multiplies(lhs, rhs):
  return [l * r for l, r in zip(lhs, rhs)]

def list_divides(lhs, rhs):
  return [l / r for l, r in zip(lhs, rhs)]

def list_modulus(lhs, rhs):
  return [l % r for l, r in zip(lhs, rhs)]

def list_dot_product(lhs, rhs=0):
  if (rhs == 0): rhs = lhs
  result = 0
  for l, r in zip(lhs, rhs): result += l * r
  return result

def str_ev(EV):
  return "[%d,%d,%d]" % EV

###def fract_2_dec(fraction):
###  list = fraction.split('/')
###  if len(list) == 2 and list[1] != 0:
###    decimal = string.atof(list[0])/string.atof(list[1])
###  else:
###    decimal = string.atof(fraction)
###  return decimal

def rlc_RTMxAnalysis(M):
  r_info = sgtbx.rot_mx_info(M.r())
  t_info = sgtbx.translation_part_info(M)
  t_intrinsic = t_info.intrinsic_part().mod_positive().as_double()
  t_shift = t_info.origin_shift().mod_positive().as_double()

  #End = list_plus(Start + map(None,r_info.ev()))
####debug
###  trans = 0
###  length = 0
####debug

  #if (r_info.type() == 1):
  if (r_info.type() < 2):
    #(rt, start, end) = ('1',(0,0,0),(0,0,0))
    return None
  #elif (r_info.type() == -1):
  #  (rt, start, end) = (str(r_info.type()),t_shift,())
  elif (abs(r_info.type()) == 2):
    trans = reduce(lambda x,y:x+y,t_intrinsic)
    if trans == 0:
      maxr = max([abs(x) for x in r_info.ev()])
      r = [float(x)/maxr for x in r_info.ev()]
      (rt, start, end) = (str(r_info.type()),t_shift,tuple(list_plus(t_shift,r)))
      #(rt, start, end) = (str(r_info.type()),t_shift,tuple(list_plus(t_shift,r_info.ev())))
    else:
      maxr = max([abs(x) for x in r_info.ev()])
      r = [float(x)/maxr for x in r_info.ev()]
      (rt, start, end) = (str(r_info.type())+"^1",t_shift,tuple(list_plus(t_shift,r)))
      #(rt, start, end) = (str(r_info.type())+"^1",t_shift,tuple(list_plus(t_shift,r_info.ev())))
  elif (r_info.type() == 3):
    if (r_info.sense() >= 0) :
      # ignore opposite sense of rotation axes since they superimpose
      trans = N.sqrt(reduce(lambda x,y:x+y,(map(lambda x,y:(y-x)*(y-x),(0,0,0),t_intrinsic))))
#      trans = N.sqrt(t_intrinsic[0]**2 + t_intrinsic[1]**2 + t_intrinsic[2]**2)
      if trans == 0:
        maxr = max([abs(x) for x in r_info.ev()])
        r = [float(x)/maxr for x in r_info.ev()]
# fudge to make sure that PyMOL actually draws the axis (move it slightly off [1,-1,1]) !!!
        r[0] = r[0]*1.000001
        (rt, start, end) = (str(r_info.type()),t_shift,tuple(list_plus(t_shift,r)))
        #(rt, start, end) = (str(r_info.type()),t_shift, tuple(list_plus(t_shift,r_info.ev())))
      else:
        maxr = max([abs(x) for x in r_info.ev()])
        r = [float(x)/maxr for x in r_info.ev()]
        #(rt, start, end) = (str(r_info.type())+ "^" + subscript ,t_shift,tuple(list_plus(t_shift,r)))
        (start, end) = (t_shift,tuple(list_plus(t_shift,r)))
        length = N.sqrt(reduce(lambda x,y:x+y,(map(lambda x,y:(y-x)*(y-x),start, end))))

#  r_info.sense() for 3^1 and 3^2 seems always to be "1" ???
#        if r_info.sense() < 0:
#          subscript = str(1-r_info.sense())
#        else:
#          subscript = str(r_info.sense())

# use ratio of trans to length to get the correct axis symbol:
# fudged the value to get the right numbers. (using length/2., rather than length/3.)
        if trans < length*0.5 :
          subscript = '1'
        else:
          subscript = '2'

        rt = str(r_info.type())+ "^" + subscript
        #(rt, start, end) = (str(r_info.type()) + "^" + subscript,t_shift, tuple(list_plus(t_shift,r_info.ev())))
###        print "Type, sense, Start, End, length, trans", rt, r_info.sense(), start, end, length, trans
#        print "type: %s, sense: %s, trans: %s, length: %s," % (r_info.type(), r_info.sense(), trans, length)
#        print "(rt, start, end)", (rt,start,end)
    else:
      return None
  #return (r_info.type(),r_info.ev(), t_intrinsic, t_shift)
  elif (r_info.sense() > 0):
    # ignore opposite sense of rotation axes since they superimpose
    trans = reduce(lambda x,y:x+y,t_intrinsic)
    if trans == 0:
      maxr = max([abs(x) for x in r_info.ev()])
      r = [float(x)/maxr for x in r_info.ev()]
      (rt, start, end) = (str(r_info.type()),t_shift,tuple(list_plus(t_shift,r)))
      #(rt, start, end) = (str(r_info.type()),t_shift, tuple(list_plus(t_shift,r_info.ev())))
    else:
      maxr = max([abs(x) for x in r_info.ev()])
      r = [float(x)/maxr for x in r_info.ev()]
      subscript =  str(int(trans*r_info.type()+.5))  # add 0.5 to fix rounding errors
      (rt, start, end) = (str(r_info.type())+ "^" + subscript ,t_shift,tuple(list_plus(t_shift,r)))
      #(rt, start, end) = (str(r_info.type()) + "^" + subscript,t_shift, tuple(list_plus(t_shift,r_info.ev())))
  #return (r_info.type(),r_info.ev(), t_intrinsic, t_shift)
  else:
    return None
#  print "type: %s, sense: %s, trans: %s, length: %s," % (r_info.type(), r_info.sense(), trans, length),
#  print "(rt, start, end)", (rt,start,end)
  return (rt, start, end)

def get_all_axes(space_group_symbol=None, space_group_info=None, extension=0):
  assert space_group_symbol is None or space_group_info is None
  shift_range = 1 # RWGK Works for the 230 reference settings; it is not
          # RWGK clear to me (rwgk) what value is needed in general.
  if (space_group_symbol is not None):
    try:
      space_group_info = sgtbx.space_group_info(symbol=space_group_symbol)
    except RuntimeError, err:
      print err
      return None

  #space_group_info.show_summary()

  axes_dict = {}
  for smx in space_group_info.group().all_ops():
    r = smx.r()
    t = smx.t()
    shift = [0,0,0]
    for shift[0] in range(-shift_range,shift_range+1):
      for shift[1] in range(-shift_range,shift_range+1):
        for shift[2] in range(-shift_range,shift_range+1):
          ts = t.plus(sgtbx.tr_vec(shift, 1)).new_denominator(t.den())
          m = sgtbx.rt_mx(r, ts)
          #print m
          rtmxanal = rlc_RTMxAnalysis(m)
          #print r, t, shift, ts, m
          if rtmxanal:
            #print rtmxanal
            axes_dict[rtmxanal] = 0
  axes_list = axes_dict.keys()
  axes_list.sort()

  # reject nonenantiomorphic space groups
  if len(axes_list) > 0 and not re.compile("[A-z]").search(space_group_symbol[1:]):
    try:
      sgtbx.space_group_info(space_group_symbol).show_summary(),
      #print len(axes_list), space_group_symbol
    except:
      print space_group, space_group_symbol
      print
      sys.exit(1)
    axes = []
    for a in axes_list:
      if len(a) == 3 and len(a[1]) == 3 and len(a[2]) == 3:
        tmp_dict = {}
        print "%4s %7.4f %7.4f %7.4f    %7.4f %7.4f %7.4f " % (a[0],a[1][0],a[1][1],a[1][2],a[2][0],a[2][1],a[2][2])
        tmp_dict['symb'] = a[0]
        start_array = N.asarray(a[1])
        end_array = N.asarray(a[2])
        start_vec = start_array - (end_array - start_array)*extension
        end_vec = end_array + (end_array - start_array)*extension
        tmp_dict['start'] = start_vec
        tmp_dict['end'] = end_vec
#rlc#        tmp_dict['start'] = a[1]
#rlc#        tmp_dict['end'] = a[2]
        axes.append(tmp_dict)
      else:
        print a
  else:
    return None

  return axes

if (__name__ == "__main__"):
  import sys
  if (len(sys.argv) == 1):
    for i in range(230):
      get_all_axes(i + 1)
  else:
    for symbol in sys.argv[1:]:
      get_all_axes(symbol)
