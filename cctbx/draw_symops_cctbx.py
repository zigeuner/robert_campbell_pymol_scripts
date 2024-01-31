#! /usr/bin/env python
# Copyright (c) 2004 Robert L. Campbell

from cctbx import uctbx, sgtbx
#import string, math
from pymol.cgo import *
from pymol import cmd

from all_axes_new import get_all_axes

import numpy as N

print "Finished importing draw_symops_cctbx.py modules"

def set_to_zero(a):
  if abs(a) < 1e-10:
    a=0
  return a

def draw_symbol(start,end,symb,color,radius=0.2):
  degtorad = N.pi/180.
  costhirty = N.cos(30.0*degtorad)
  sinthirty = N.sin(30.0*degtorad)
  cossixty = N.cos(60.0*degtorad)
  sinsixty = N.sin(60.0*degtorad)
  symb_obj = []

  if symb == '2' or symb == '2^1':
    pass

  elif symb == '3' or symb == '3^1' or symb == '3^2':
    symb_obj = [ BEGIN, TRIANGLES, COLOR ] + color
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([start]) + N.array([radius, 0, 0]))[0].tolist()
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([start]) + N.array([-radius*sinthirty, radius*costhirty, 0]))[0].tolist()
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([start]) + N.array([-radius*sinthirty, -radius*costhirty, 0]))[0].tolist()

    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([end]) + N.array([radius, 0, 0]))[0].tolist()
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([end]) + N.array([-radius*sinthirty, radius*costhirty, 0]))[0].tolist()
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([end]) + N.array([-radius*sinthirty, -radius*costhirty, 0]))[0].tolist()
    symb_obj.append(END)

  elif symb == '4' or symb == '4^1' or symb == '4^2' or symb == '4^3':
    symb_obj = [ BEGIN, TRIANGLES, COLOR ] + color
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([start]) + N.array([radius, radius, 0]))[0].tolist()
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([start]) + N.array([-radius, radius, 0]))[0].tolist()
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([start]) + N.array([-radius, -radius, 0]))[0].tolist()

    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([start]) + N.array([radius, radius, 0]))[0].tolist()
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([start]) + N.array([-radius, -radius, 0]))[0].tolist()
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([start]) + N.array([radius, -radius, 0]))[0].tolist()

    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([end]) + N.array([radius, radius, 0]))[0].tolist()
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([end]) + N.array([-radius, radius, 0]))[0].tolist()
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([end]) + N.array([-radius, -radius, 0]))[0].tolist()

    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([end]) + N.array([radius, radius, 0]))[0].tolist()
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([end]) + N.array([-radius, -radius, 0]))[0].tolist()
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([end]) + N.array([radius, -radius, 0]))[0].tolist()
    symb_obj.append(END)

  elif symb == '6' or symb == '6^1' or symb == '6^2' or symb == '6^3' or symb == '6^4' or symb == '6^5':
    symb_obj = [ BEGIN, TRIANGLES, COLOR ] + color
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([start]) + N.array([radius, 0, 0]))[0].tolist()
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([start]) + N.array([radius*cossixty, radius*sinsixty, 0]))[0].tolist()
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([start]) + N.array([-radius*cossixty, radius*sinsixty, 0]))[0].tolist()

    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([start]) + N.array([radius, 0, 0]))[0].tolist()
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([start]) + N.array([-radius*cossixty, radius*sinsixty, 0]))[0].tolist()
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([start]) + N.array([-radius, 0, 0]))[0].tolist()

    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([start]) + N.array([radius, 0, 0]))[0].tolist()
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([start]) + N.array([-radius, 0, 0]))[0].tolist()
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([start]) + N.array([-radius*cossixty, -radius*sinsixty, 0]))[0].tolist()

    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([start]) + N.array([radius, 0, 0]))[0].tolist()
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([start]) + N.array([-radius*cossixty, -radius*sinsixty, 0]))[0].tolist()
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([start]) + N.array([radius*cossixty, -radius*sinsixty, 0]))[0].tolist()

    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([end]) + N.array([radius, 0, 0]))[0].tolist()
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([end]) + N.array([radius*cossixty, radius*sinsixty, 0]))[0].tolist()
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([end]) + N.array([-radius*cossixty, radius*sinsixty, 0]))[0].tolist()

    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([end]) + N.array([radius, 0, 0]))[0].tolist()
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([end]) + N.array([-radius*cossixty, radius*sinsixty, 0]))[0].tolist()
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([end]) + N.array([-radius, 0, 0]))[0].tolist()

    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([end]) + N.array([radius, 0, 0]))[0].tolist()
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([end]) + N.array([-radius, 0, 0]))[0].tolist()
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([end]) + N.array([-radius*cossixty, -radius*sinsixty, 0]))[0].tolist()

    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([end]) + N.array([radius, 0, 0]))[0].tolist()
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([end]) + N.array([-radius*cossixty, -radius*sinsixty, 0]))[0].tolist()
    symb_obj.append(VERTEX)
    symb_obj = symb_obj + (N.array([end]) + N.array([radius*cossixty, -radius*sinsixty, 0]))[0].tolist()

    symb_obj.append(END)

  return symb_obj

def draw_symops(obj,radius=0.2,extension=0,prefix='',individual_axes=0):
  """
  From pymol issue the "run draw_symops_cctbx.py" command to load the script,
  then issue the "draw_symops(object,<optional radius>,<optional extension>)" command
  to actually run it and create the cgo object.

  e.g. load 1avv.pdb
       run draw_symops_cctbx.py
       draw_symops 1avv, 0.5, .2
         or draw_symops('1avv',.5,.2)
         or draw_symops 1avv, radius=.5, extension=.2
         or draw_symops 1avv, radius=.5, extension=.2, prefix='1avv_'

  The different axis types appear as different objects on the PyMOL menu so they can be turned
  on and off individually.

  The 'extension' parameter is a fractional increase in the length of
  each symmetry operator axis drawn. i.e. a value of 0 is the default
  and a value of .2 increases the length by 20% at each end.

  The 'prefix' parameter allows you to include your own prefix into the
  name of the space group operator objects.

  The 'individual_axes' option selects whether to use a separate object
  for each axis (individual_axes=1), or whether to group them by
  symmetry symbol (default, individual_axes=0).

  See also help(draw_symops_param) to draw operators by specifying the
  space group and cell dimensions directly (i.e. not loaded from a pdb
  file).  The same options apply to that function.

  """
  radius=float(radius)
  extension=float(extension)
  cell_info=cmd.get_symmetry(obj)
  draw_symops_param(cell_info[0:6],cell_info[6],radius,extension,prefix,individual_axes)
  individual_axes=int(individual_axes)

def draw_symops_param(cell_param_list,sg,radius=0.2,extension=0,prefix='',individual_axes=0):
  """
  If you wish to draw the symmetry operators for any cell without the need to load a
  pdb file, then do this:

  e.g. run draw_symops_cctbx.py
       draw_symops_param((45.2,45.2,70.8,90.,90.,120.),'p3121',0.5,0.1)

  to generate the symmetry operators for this trigonal space group "p 31 2 1"
  of radius .5 with 10% added as an extension at each end.
  """
  radius=float(radius)
  extension=float(extension)
  individual_axes=int(individual_axes)

  U=uctbx.unit_cell((cell_param_list))

  """
#rotation axes
#    "2" "yellow",
#    "3" "orange",
#    "4" "mauve",
#    "6" "purple",

#screw axes (all sub_1 axes are green)
#    "21" "green",
#    "31" "green",
#    "32" "lime",
#    "41" "green",
#    "42" "cyan",
#    "43" "iceblue",
#    "61" "green",
#    "62" "silver",
#    "63" "cyan",
#    "64" "iceblue",
#    "65" "blue",
  """
  color = {
    "2" : [1.0, 1.0, 0.0],
    "3" : [1.0, 0.5, 0.0],
    "4" : [1.0, 0.5, 1.0],
    "6" : [1.0, 0.0, 1.0],
    "2^1" : [0.0, 1.0, 0.0],
    "3^1" : [0.0, 1.0, 0.0],
    "3^2" : [0.5, 1.0, 0.5],
    "4^1" : [0.0, 1.0, 0.0],
    "4^2" : [0.0, 1.0, 1.0],
    "4^3" : [0.5, 0.5, 1.0],
    "6^1" : [0.0, 1.0, 0.0],
    "6^2" : [0.8, 0.8, 0.8],
    "6^3" : [0.0, 1.0, 1.0],
    "6^4" : [0.5, 0.5, 1.0],
    "6^5" : [0.0, 0.0, 1.0],
    }

  sg = sg.upper()
  symop_axes = get_all_axes(sg,extension=extension)

  #CYLINDER = 'CYLINDER'
  ax_obj = {}
  ax_list_obj = []

  #debug_out = open('debug.log','w')

  if symop_axes:
    for i in range(len(symop_axes)):
      #print symop_axes[i]
      start = map(set_to_zero,U.orthogonalize(map(None,symop_axes[i]['start'])))
      end = map(set_to_zero,U.orthogonalize(map(None,symop_axes[i]['end'])))
###############################################################################
# Tried rounding off start and end values in order to understand why axes go
# missing in the drawing, but seem to be present in the cgo.  Doesn't help!
# e.g. for space group 'p23' one of the 3-fold rotations is missing (0,0,0 -> x,-x,x)
# changing one cell axis to something ever so slightly different recovers the axis
# e.g. set cell to be (30.00001,30.,30.,90.,90.,90) and it works!
#    start = map(lambda x: round(x,3),U.orthogonalize(symop_axes[i]['start']))
#    end = map(lambda x: round(x,3),U.orthogonalize(symop_axes[i]['end']))
###############################################################################
      color_ax = color[symop_axes[i]['symb']]
      symb_ax = symop_axes[i]['symb']

      #print "axis: ",symb_ax, start, end
# select whether to have an individual object for each axis or whether to group them
# by symb_ax name
      if individual_axes:
        symb_ax_key = str(i) + '_' + symb_ax
      else:
        symb_ax_key = symb_ax

      if ax_obj.has_key(symb_ax):
        ax_obj[symb_ax_key].append(CYLINDER)
      else:
        ax_obj[symb_ax_key] = [CYLINDER]

      ax_obj[symb_ax_key] = ax_obj[symb_ax_key] + start + end + [radius]
      ax_obj[symb_ax_key] = ax_obj[symb_ax_key] + color[symb_ax] + color[symb_ax]
      ax_obj[symb_ax_key] = ax_obj[symb_ax_key] + draw_symbol(start,end,symb_ax,color[symb_ax],radius*6.)

#    #######################################################################################
#      # Debugging output to try to understand why some axes go missing in the drawing.
#      # They don't appear to be missing from the cgo object, though!
#      for xxx in ax_obj[symb_ax]:
#        if xxx == 9.0:
#          #print "\n\n",xxx
#          xxx = "\n\n" + str(xxx) + " "
#          debug_out.write(xxx)
#        else:
#          #print xxx
#          #xxx = "\n" + str(xxx) + " "
#          xxx = str(xxx) + " "
#          debug_out.write(xxx)
#        #print ax_obj[symb_ax]
#    debug_out.write("\n\n")
#    big_string = str(ax_obj)
#    debug_out.write(big_string)
#    # End of debugging output
#    #######################################################################################

  else:
    print "\nNo symmetry axes found for this space group: %s\n" % sg

  key_list = ax_obj.keys()
# if drawing individual objects sort key_list based on counter (what comes before first '_')
  if individual_axes:
    key_list.sort(lambda x,y: (cmp(int(x.split('_')[0]),int(y.split('_')[0]))))
# else do a plain text sort
  else:
    key_list.sort()

  for key in key_list:
    name = prefix + sg + "_" + str(key)
    cmd.load_cgo(ax_obj[key],name)
    #debug_out.write("\n\n" + key + "\n" + str(ax_obj[key]))
  #return ax_obj

cmd.extend("draw_symops",draw_symops)
cmd.extend("draw_symops_param",draw_symops_param)
