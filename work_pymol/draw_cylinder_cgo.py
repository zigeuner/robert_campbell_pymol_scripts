#! /usr/bin/python
# Copyright (c) 2010 Robert L. Campbell
#
# A PyMOL script for drawing a CGO cylinder from the coordinates of three atoms (pk1 and pk2 by default)

from pymol.cgo import *
from pymol import cmd


def draw_cylinder_cgo(name,start,end,radius=0.1,start_color=[1,1,0],end_color=[1,1,0]):

  """
DESCRIPTION
    Draw a CGO cylinder between two arbitrary coordinates

USAGE
    draw_cylinder_cgo name, start, end, radius, start_color, end_color

    where start and end are 3-element vectors and the colors are either a name or
    a 3-element RGB list enclosed in parenthesis defining the color of the cylinder
    end (where each value of R, G and B is between 0 and 1 inclusive).

    The radius defaults to 0.1A and the colors to yellow: (1,1,0)

  """

  # Convert args to floating point numbers
  start = list(map(float,start))
  end = list(map(float,end))
  radius = float(radius)
  if type(start_color) is str:
    try:
      tup_start_color = list(map(float, start_color.replace('(','').replace(')','').split(',')))
    except ValueError:
#    print "converting start_color %s to list" % (start_color)
      tup_start_color = list(cmd.get_color_tuple(start_color))
#    print "Conversion for ", start_color," to ", tup_start_color
  elif type(start_color) is list or type(start_color) is tuple:
    tup_start_color = list(start_color)

  if type(end_color) is str:
    try:
      tup_end_color = list(map(float, end_color.replace('(','').replace(')','').split(',')))
    except ValueError:
#    print "converting end_color %s to list" % (end_color)
      tup_end_color = list(cmd.get_color_tuple(end_color))
#    print "Conversion for ", end_color," to ", tup_end_color
  elif type(end_color) is list or type(end_color) is tuple:
    tup_end_color = list(end_color)

  # Create the CGO objects
  obj = [CYLINDER]
  obj +=  start + end + [radius] + tup_start_color + tup_end_color

  # Display it
  cmd.load_cgo(obj,name)

def draw_cylinder(name,atom1='(pk1)',atom2='(pk2)',radius=0.1,start_color=[1,1,0],end_color=[1,1,0]):
  """
DESCRIPTION
    Draw a CGO cylinder between two atoms

USAGE
    draw_cylinder name, atom1, atom2, radius, start_color, end_color

    where each atom is a standard PyMOL selection (defaults to pk1 and pk2)
    and color is a name or a 3-element RGB tuple defining the color
    of the cylinder (where each value of R, G and B is between 0 and 1
    inclusive).  The radius defaults to 0.1A and the colors default
    to yellow (1,1,0).

  """

  coor1 = cmd.get_model(atom1).atom[0].coord
  coor2 = cmd.get_model(atom2).atom[0].coord
  draw_cylinder_cgo(name,coor1,coor2,radius,start_color,end_color)

cmd.extend("draw_cylinder", draw_cylinder)
cmd.extend("draw_cylinder_cgo", draw_cylinder_cgo)
