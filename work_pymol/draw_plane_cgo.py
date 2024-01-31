#! /usr/bin/python
# Copyright (c) 2010 Robert L. Campbell
#
# A PyMOL script for drawing a CGO plane from the coordinates of three atoms (pk1,pk2,pk3 by default)

from pymol.cgo import *
from pymol import cmd


# A helper function for computing the normal to a triangular facet
def compute_normal(x1, y1, z1, x2, y2, z2, x3, y3, z3):

  nx = (y2-y1)*(z3-z2) - (z2-z1)*(y3-y2)
  ny = (z2-z1)*(x3-x2) - (x2-x1)*(z3-z2)
  nz = (x2-x1)*(y3-y2) - (y2-y1)*(x3-x2)

  return (nx,ny,nz)



def draw_plane_cgo(name,apex1,apex2,apex3,apex4,color=(1,1,1)):

  """
DESCRIPTION
    Create a CGO plane from three arbitary coordinates

USAGE
    draw_plane_cgo apex1, apex2, apex3, apex4, color

    where each apex is a 3-element vector and color is a 3-element RGB
    list defining the color of the plane (where each value of R, G
    and B is between 0 and 1 inclusive).

  """

  # Convert args to floating point numbers
  x1,y1,z1 = map(float,apex1)
  x2,y2,z2 = map(float,apex2)
  x3,y3,z3 = map(float,apex3)
  x4,y4,z4 = map(float,apex4)
  if type(color) == type(''):
    color = map(float,color.replace('(','').replace(')','').split(','))

  # Compute the normal vector for the triangle
  normal1 = compute_normal(x1, y1, z1, x2, y2, z2, x3, y3, z3)
  normal2 = compute_normal(x1, y1, z1, x3, y3, z3, x4, y4, z4)
  normal3 = compute_normal(x2, y2, z2, x3, y3, z3, x4, y4, z4)

  # Create the CGO objects
  obj = [

    BEGIN, TRIANGLE_STRIP,

    COLOR, color[0], color[1], color[2],
    NORMAL, normal1[0], normal1[1], normal1[2],
    VERTEX, x1, y1, z1,
    VERTEX, x2, y2, z2,
    VERTEX, x3, y3, z3,
    VERTEX, x4, y4, z4,

    END
  ]

  # Display them
  cmd.load_cgo(obj,name)

def draw_plane(name,atom1='(pk1)',atom2='(pk2)',atom3='(pk3)',atom4='(pk4)',color=(1,1,1)):
  """
DESCRIPTION
    Create a CGO plane from four atomic coordinates

USAGE
    draw_plane name, atom1, atom2, atom3, atom4, color

    where each atom is a standard PyMOL selection (defaults to pk1,pk2
    and pk3) and color is a 3-element RGB tuple defining the color
    of the plane (where each value of R, G and B is between 0 and 1
    inclusive).  The color defaults to (1,1,1).

  """

# get coordinates from atom selections
  coor1 = cmd.get_model(atom1).atom[0].coord
  coor2 = cmd.get_model(atom2).atom[0].coord
  coor3 = cmd.get_model(atom3).atom[0].coord
  coor4 = cmd.get_model(atom4).atom[0].coord
  draw_plane_cgo(name,coor1,coor2,coor3,coor4,color)

cmd.extend("draw_plane", draw_plane)
