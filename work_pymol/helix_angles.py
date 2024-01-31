#! /usr/bin/env python3
# coding=latin-1


import numpy,sys,getopt

# magnitude of a vector
def magvect(v):
  magnitude = numpy.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
  return magnitude


# angle between vectors
def angle(v1,v2):
  ang = numpy.arccos(numpy.dot(v1,v2)/(magvect(v1)*magvect(v2)))
  ang = ang*180/numpy.pi
#debug
#  print("magv1 magv2 dv1v2 cos(angle)", magv1, magv2, dv1v2, dv1v2/(magv1*magv2))
  return ang


# distance between vectors
def distance(start1,v1,start2,v2):
  cross_prod = numpy.cross(v1,v2)
  mx = magvect(cross_prod)
  norm = cross_prod/mx
  diff = start1 - start2
  dist = numpy.fabs(numpy.dot(norm,diff))
  return dist


# least-squares fit of two arrays, x and y, and return first and last points from fit
def lsq(x,y):
  slope,intcpt =  numpy.polyfit(x,y,1,full=False)
  first_point = slope*x[0]+intcpt
  last_point = slope*x[-1]+intcpt
  return first_point,last_point

def usage():
  print("""
  helix_angles.py [filename chain1 start1 end1 chain2 start2 end2]
  """)
  sys.exit(0)

def get_helices_from_PDB(input_data):
  try:
    filename,chain1,res1a,res1b,chain2,res2a,res2b = input_data
  except ValueError:
    filename,res1a,res1b,res2a,res2b = input_data
    chain1=' '
    chain2=' '

  res1a = int(res1a)
  res1b = int(res1b)
  res2a = int(res2a)
  res2b = int(res2b)

  i=0
  j=0
  count1 = []
  count2 = []
  coord1 = []
  coord2 = []
#  x1 = []
#  y1 = []
#  z1 = []
#  x2 = []
#  y2 = []
#  z2 = []

  print("Opening file: ",filename)
  print("Looking for residues %d and %d in chain %s and for residues %d and %d in chain %s." % (res1a,res1b,chain1,res2a,res2b,chain2))

  Infile = open(filename,'r')
  for line in Infile.readlines():

    if (line[0:4] == 'ATOM' and 
           (line[12:16] == ' CA ' or
            line[12:16] == ' N  ' or
            line[12:16] == ' C  ' or
            line[12:16] == ' P  ') or
       (line[0:6] == 'HETATM' and line[12:16] == ' P  ')):
  # Find C-alphas for protein and Phosphates for DNA
      chain = line[21:22]
      res = int(line[22:26])

      # Find the residues we want in the first helix
      if (chain == chain1 and (res >= res1a and res <= res1b)):
        i += 1
        count1.append(i)
#        x1.append(float(line[30:38]))
#        y1.append(float(line[38:46]))
#        z1.append(float(line[46:54]))
        coord1.append((float(line[30:38]),float(line[38:46]),float(line[46:54])))
# for debugging
#        print("i,chain,res,count,xyz ", i,chain,res, count1[i],x1[i],y1[i],z1[i])


      # Find the residues we want in the second helix
      elif (chain == chain2 and (res >= res2a and res <= res2b)):
        j += 1
        count2.append(j)
#        x2.append(float(line[30:38]))
#        y2.append(float(line[38:46]))
#        z2.append(float(line[46:54]))
        coord2.append((float(line[30:38]),float(line[38:46]),float(line[46:54])))
# for debugging
#        print("j,chain,res,count,xyz ", j,chain,res, count2[j],x2[j],y2[j],z2[j])

  nres1 = i
  nres2 = j
  print("Found", nres1, "and", nres2, "atoms in the two helices, respectively.")
#debug
#  for i in range(nres1):
#    print(count1[i],x1[i],y1[i],z1[i])
#  for i in range(nres2):
#    print(count2[i],x2[i],y2[i],z2[i])
  Infile.close()
  return coord1,coord2


def helix_angle(selection1,selection2,color=None,radius=0.2,object_name=None):
  """
  AUTHOR
    Robert L. Campbell

    Calculate angles between helix axes (and optionally draw axes as cylinders)

  USAGE
    helix_angle(selection1='sel',selection2='sel',color='red',radius=0.2,object_name='helix_axes')

    - color defaults to red
    - color can be specified by name, or by color tuple with values
      between 0 and 1, e.g. red is (1,0,0) and magenta is (1,0,1)
    - radius defaults to 0.2A
    - object_name defaults to 'helix_axes'.
  """
  from pymol import cmd
  from pymol import cgo

# set defaults

  if not color:
    tup_color = [1.,0.,0.]
  if type(color) is str:
    try:
      tup_color = list(map(float, color.replace('(','').replace(')','').split(',')))
    except ValueError:
#    print("converting color %s to list" % (color))
      tup_color = list(cmd.get_color_tuple(color))
#    print("Conversion for ", color," to ", tup_color)
  elif type(color) is list or type(color) is tuple:
    tup_color = list(color)

    radius = float(radius)

# get atoms in the two (helical) selections
  print("Looking for residues in %s and %s." % (selection1,selection2))
  m1 = cmd.get_model(selection1)
  m2 = cmd.get_model(selection2)
  coord1 = []
  coord2 = []

  if len(m1.atom) == 0:
    print("Sorry, no atoms selected in selection 1: ", selection1)
    sys.exit(1)
  elif len(m2.atom) == 0:
    print("Sorry, no atoms selected in selection 2: ", selection2)
    sys.exit(1)

  else:
    if len(m1.atom) > 1 or len(m2.atom) > 1:
      for a in m1.atom:
        if (a.name == 'CA' or a.name == 'N' or a.name == 'C' or a.name == 'P'):
          coord1.append(a.coord)

      for b in m2.atom:
        if (b.name == 'CA' or b.name == 'N' or b.name == 'C' or b.name == 'P'):
          coord2.append(b.coord)
    else:
      print("Sorry, there must be more than one atom in each selection in order to define the axes")
      sys.exit(1)

  print("Found", len(coord1), "and", len(coord2), "atoms in the two helices, respectively.")

  s1,v1,s2,v2 = calc_helix_axes((coord1,coord2))
  if object_name:
    from pymol import cgo
    cyl_obj = draw_vectors(((s1,v1),(s2,v2)),radius,color=tup_color)
    cmd.load_cgo(cyl_obj,object_name)

  return coord1,coord2

def calc_helix_axes(xxx_todo_changeme):
  # Do lsq fitting of x, y and z values for each helix
  (coord1,coord2) = xxx_todo_changeme
  if (len(coord1) > 0 and len(coord2) > 0):
    coord1 = numpy.asarray(coord1)
    coord2 = numpy.asarray(coord2)
    count1 = len(coord1)
    count2 = len(coord2)

    start1 = []
    end1 = []
    vect1 = []
    start2 = []
    end2 = []
    vect2 = []

    (s,e) = lsq(list(range(count1)),coord1[:,0])
    start1.append(s)
    end1.append(e)
    (s,e) = lsq(list(range(count1)),coord1[:,1])
    start1.append(s)
    end1.append(e)
    (s,e) = lsq(list(range(count1)),coord1[:,2])
    start1.append(s)
    end1.append(e)

    vect1.append(end1[0] - start1[0])
    vect1.append(end1[1] - start1[1])
    vect1.append(end1[2] - start1[2])

    (s,e) = lsq(list(range(count2)),coord2[:,0])
    start2.append(s)
    end2.append(e)
    (s,e) = lsq(list(range(count2)),coord2[:,1])
    start2.append(s)
    end2.append(e)
    (s,e) = lsq(list(range(count2)),coord2[:,2])
    start2.append(s)
    end2.append(e)

    vect2.append(end2[0] - start2[0])
    vect2.append(end2[1] - start2[1])
    vect2.append(end2[2] - start2[2])

    start1 = numpy.asarray(start1)
    start2 = numpy.asarray(start2)
    vect1 = numpy.asarray(vect1)
    vect2 = numpy.asarray(vect2)

    print("\nStart 1 = (%9.4f, %9.4f, %9.4f)" % (start1[0], start1[1], start1[2]), \
    "Vector 1 = (%9.4f, %9.4f, %9.4f)" % (vect1[0], vect1[1], vect1[2]))

    print("Start 2 = (%9.4f, %9.4f, %9.4f)" % (start2[0], start2[1], start2[2]), \
    "Vector 2 = (%9.4f, %9.4f, %9.4f)" % (vect2[0], vect2[1], vect2[2]))

    ang1_2 = angle(vect1,vect2)
    dist1_2 = distance(start1,vect1,start2,vect2)

    print("\nAngle between helix axes: %6.2f degrees \nDistance between helix axes: %6.2f A\n" % (ang1_2,dist1_2))

    return start1,vect1,start2,vect2

  else:
    print("Missing coordinates for one or both helices.  Check your input.")

def draw_vectors(vectors,radius,color=(1.0,0.0,0.0)):
  """
  vectors is an array of start and end pairs of vectors
  """
  radius = float(radius)
  size = radius*5.
  cone_height = radius*5
  cone_base_radius = radius*2
  obj = []

  for s,e in vectors:
#    CYLINDER='CYLINDER'
#    CONE='CONE'
#    CYLINDER=9.0
#    CONE=27.0
    CYLINDER=cgo.CYLINDER
    CONE=cgo.CONE
    norm = numpy.linalg.norm(numpy.subtract(e,s))
    obj.append(CYLINDER)
    obj.extend(list(s))
    obj.extend(list(s+e))
    obj.append(radius)
    obj.extend(color)
    obj.extend(color)
#    obj.append(CONE)
#    obj.extend(list(e))
#    obj.extend(numpy.add(e,norm*cone_height))
#    obj.append(cone_base_radius)
#    obj.append(0.0)
#    obj.extend(color)
#    obj.extend([1.0,1.0])
#  print(obj)
  return obj

if __name__ == "__main__":

  #command line stuff here

  create_cgo_object = 0
  try:
    opts,args = getopt.getopt(sys.argv[1:],'hc',['help','cgo'])
  except:
    sys.stderr.write("\n*********************************\n")
    sys.stderr.write("\n      Unknown options %s\n" % str(sys.argv[1:]))
    sys.stderr.write("\n*********************************\n\n")
    usage()

  for o,a in opts:
    if o in ('-h', '--help'):
      usage()
    if o in ('-c', '--cgo'):
      create_cgo_object = 1

  if len(args) == 0:
      #    file_in = args[0]

    print("Enter filename, and starting and ending residue numbers, with ")
    print("appropriate chain identifier, for both helices.")
    print("<E.g.  prot_dna.pdb A 15 30 D 1 12>\n")
    for data in sys.stdin.readlines():
      s1,v1,s2,v2 = calc_helix_axes(get_helices_from_PDB(data.split()))

      if create_cgo_object:
        draw_vectors(((s1,v1),(s2,v2)),radius=0.5,color=(1.0,0.,0.,1.0,0.,0.))

  else:

    s1,v1,s2,v2 = calc_helix_axes(get_helices_from_PDB(args))

    if create_cgo_object:
      draw_vectors(((s1,v1),(s2,v2)),radius=0.3,color=(1.0,0.,0.,1.0,0.,0.))

else:
  from pymol import cmd
  cmd.extend("helix_angle",helix_angle)

# End of functions
