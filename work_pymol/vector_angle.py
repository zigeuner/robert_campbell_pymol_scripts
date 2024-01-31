#! /usr/bin/env python
# coding=latin-1 
# Copyright (c) 2008 Robert L. Campbell (rlc1@queensu.ca)

import numpy,sys

# magnitude of a vector
def magvect(v):
  magnitude = numpy.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
  return magnitude


# angle between vectors
def angle(v1,v2):
  ang = numpy.arccos(numpy.dot(v1,v2)/(magvect(v1)*magvect(v2)))
  ang = ang*180/numpy.pi
#  if (ang > 90):
#    ang = 180 - ang
    #    print("magv1 magv2 dv1v2 cos(angle)", magv1, magv2, dv1v2, dv1v2/(magv1*magv2))
  return ang


# distance between vectors
def distance(start1,v1,start2,v2):
  cross_prod = numpy.cross(v1,v2)
  mx = magvect(cross_prod)
  norm = cross_prod/mx
  diff = start1 - start2
  dist = numpy.fabs(numpy.dot(norm,diff))
  return dist

def get_vector_angle(start1,end1,start2,end2):
  """
  get_vector_angle start1, end1, start2, end2
  where:
    start1,end1 = atom selections for first vector
    start2,end2 = atom selections for second vector
  """
  #from pymol import cmd

# get atoms from models for each selection
  s1 = cmd.get_model(start1)
  name_start1 = cmd.get_object_list(start1)[0]
  e1 = cmd.get_model(end1)
  name_end1 = cmd.get_object_list(end1)[0]
  s2 = cmd.get_model(start2)
  name_start2 = cmd.get_object_list(start2)[0]
  e2 = cmd.get_model(end2)
  name_end2 = cmd.get_object_list(end2)[0]

  atoms = []
  for mol in [s1,e1,s2,e2]:
    if len(mol.atom) > 1:
      print("You must select only one atom per selection")
    else:
      a = mol.atom[0]
      atoms.append((a.chain,a.resn,a.resi,a.coord))

  vector1 = numpy.subtract(atoms[1][3],atoms[0][3])
  start1 = numpy.asarray(atoms[0][3])
  end1 = numpy.asarray(atoms[1][3])
  vector2 = numpy.subtract(atoms[3][3],atoms[2][3])
  start2 = numpy.asarray(atoms[2][3])
  end2 = numpy.asarray(atoms[3][3])

  print("Calculating angle and distance between the vectors from:")
  print(" ", name_start1,s1.atom[0].chain, s1.atom[0].resn,s1.atom[0].resi,s1.atom[0].name, "to ",end=" ")
  print(name_end1,e1.atom[0].chain, e1.atom[0].resn,e1.atom[0].resi,e1.atom[0].name)
  print("and")
  print(" ", name_start2,s2.atom[0].chain, s2.atom[0].resn,s2.atom[0].resi,s2.atom[0].name, "to ",end=" ")
  print(name_end2,e2.atom[0].chain, e2.atom[0].resn,e2.atom[0].resi,e2.atom[0].name)
  print("")
  print("Angle between    [%.3f %.3f %.3f]" % (vector1[0],vector1[1],vector1[2]),end=" ")
  print("and [%.3f %.3f %.3f] " % (vector2[0],vector2[1],vector2[2]), " is %.2f degrees" % angle(vector1,vector2))
  print("Distance between [%.3f %.3f %.3f]" %  (vector1[0],vector1[1],vector1[2]),end=" ")
  print("and [%.3f %.3f %.3f] " % (vector2[0],vector2[1],vector2[2]), " is %.2f " % distance(start1,vector1,start2,vector2))



if __name__ == "__main__":

  points = []
  vector = []
  start = []
  end = []
  for input in sys.stdin.readlines():
    if len(input.split()) == 6:
      points = list(map(float,input.replace(',',' ').split()))
      vector.append(numpy.subtract(points[3:6],points[0:3]))
      start.append(numpy.asarray(points[0:3]))
      end.append(numpy.asarray(points[3:6]))
    elif len(input.replace(',',' ').split()) == 3:
      vector.append(list(map(float,input.split)))
    else:
      print("Enter either 2 3D coordinates (i.e. 6 values) or 1 3D vector (3 values).")
      sys.exit(1)

  for i in range(0,len(vector),2):
      print("Angle between ", vector[i], " and ",vector[i+1], " is %.2f degrees" % angle(vector[i],vector[i+1]))
      if len(start) > 0:
        print("Distance between ", vector[i], " and ",vector[i+1], " is %.2f " % distance(start[i],vector[i],start[i+1],vector[i+1]))

else:
  from pymol import cmd
  cmd.extend("get_vector_angle",get_vector_angle)
