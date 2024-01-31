#! /usr/bin/python
# rlc distance_states version 1.0
# Copyright (c) 2006 Robert L. Campbell

# don't bother with this, just include the distance function right here
#import MyPDB
from pymol import cmd
import math

def dist(x1,x2):
  """
  Calculate the distance between two coordinates.
  Returns a float
  """
  return math.sqrt((x1[0]-x2[0])**2 + (x1[1]-x2[1])**2 + (x1[2]-x2[2])**2)
            

def distance_states(selection1,selection2,same=0):
  """
  AUTHOR 
      
    Robert L. Campbell 


  USAGE

    distance_states selection1,selection2,[same=1]

    Each selection must be only one atom, but it is assumed (same=0)
    that one or both belong to different objects with multiple states.
    The function will loop over all states of each object and print out
    the distances between the pair of atoms specified.

    One can simply select the atoms (as pk1 and pk2) and then specify that on the command line:

      distance_states pk1,pk2

    or one can be explicit:

      distance_states object1 and resi 25 and name OH, object2 and resi 13 and name N

   If you want to measure a distance within an object for each state of that object,
   of between two object that have the same number of states, then pick the atoms and do:

     distance_states pk1,pk2,same=1

  """
# make sure "same" flag is integer
  same = int(same)

  print "Measuring distance between: ",selection1
  print "                       and: ",selection2

  print "State1 State2   Distance"
# states start counting at 1, so need to add 1 to the range
  for m in range(1,cmd.count_states(selection1)+1):
    sel1 = cmd.get_model(selection1,state=m)
    if len(sel1.atom) > 1:
      print "Must select only one atom for selection1"
    else:
      coord1 = sel1.atom[0].coord

# if same is 0, false, then both selections do not belong to the same object
# and therefore we need to loop over the states in selection2
    if same == 0:
      for n in range(1,cmd.count_states(selection2)+1):
        sel2 = cmd.get_model(selection2,state=n)
        if len(sel2.atom) > 1:
          print "Must select only one atom for selection2"
        else:
          coord2 = sel2.atom[0].coord
          print "%5d %5d     %6.3f" % (m,n,dist(coord1,coord2))
# same is not 0, so we assume that we want to measure the distance between
# selection1 and selection2 in the same state number only
    else:
      sel2 = cmd.get_model(selection2,state=m)
      if len(sel2.atom) > 1:
        print "Must select only one atom for selection2"
      else:
        coord2 = sel2.atom[0].coord
        print "%5d %5d     %6.3f" % (m,m,dist(coord1,coord2))

cmd.extend("distance_states",distance_states)
