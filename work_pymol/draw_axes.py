# Copyright (c) 2003 Robert L. Campbell

from pymol.cgo import *
from pymol import cmd
from pymol.vfont import plain

# create the axes object, draw axes with cylinders coloured red, green,
#blue for X, Y and Z

def draw_axes(radius=.2,length=2.,name="axes"):
  """
  usage: draw_axes radius=.05,length=2.,name=axes
  or draw_axes, .1, 10, axes
  default radius is 0.2, default length is 2 A, default name is "axes"

  """
  radius = float(radius)
  length = float(length)
  size = radius*5.
  cone_height = radius*5
  cone_base_radius = radius*2

  origin_offset = radius * -10.
  obj = [
     CYLINDER, 0., 0., 0., length, 0., 0., radius, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
     CYLINDER, 0., 0., 0., 0., length, 0., radius, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
     CYLINDER, 0., 0., 0., 0., 0., length, radius, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
     CONE, length, 0., 0., cone_height+length, 0., 0., cone_base_radius, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0,
     CONE, 0., length, 0., 0., cone_height+length, 0., cone_base_radius, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0,
     CONE, 0., 0., length, 0., 0., cone_height+length, cone_base_radius, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0,
     ]

# add labels to axes object
  text = [COLOR, 1.0, 1.0, 0.0,]

  cyl_text(text,plain,[origin_offset,origin_offset,-1],'Origin',radius,axes=[[size,0.0,0.0],[0.0,size,0.0],[0.0,0.0,size]])
  cyl_text(text,plain,[length+(2*radius),0.,0.],'X',radius,axes=[[size,0.0,0.0],[0.0,size,0.0],[0.0,0.0,size]])
  cyl_text(text,plain,[0.,length+(2*radius),0.],'Y',radius,axes=[[size,0.0,0.0],[0.0,size,0.0],[0.0,0.0,size]])
  cyl_text(text,plain,[0.,0.,length+(2*radius)],'Z',radius,axes=[[size,0.0,0.0],[0.0,size,0.0],[0.0,0.0,size]])

# then we load it into PyMOL
  #print obj
  cmd.load_cgo(obj,name)
  cmd.load_cgo(text,name + '_labels')

cmd.extend("draw_axes",draw_axes)

