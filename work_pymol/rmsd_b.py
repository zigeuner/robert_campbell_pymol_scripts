#! /usr/bin/python3
# Copyright (c) 2009 Robert L. Campbell (rlc1@queensu.ca)
import math
from pymol import cmd,stored

def distance2(x1,x2):
  """
  Calculate the square of the distance between two coordinates.
  Returns a float
  """
  dist2 = (x1[0]-x2[0])**2 + (x1[1]-x2[1])**2 + (x1[2]-x2[2])**2
  return dist2


def rmsd_b(selection1,selection2,byres=0):
  """
  AUTHOR

    Robert L. Campbell

  USAGE

    rmsd_b selection1, selection2 [,byres=0]

    Calculate the RMS deviation for each residue in a two selections and modify the B-factors
    for selection 1 to contain these RMSD values.

    The two selections must have the same number of atoms.  The easiest way to guarantee
    that is to specify an alignment object in the selection.  E.g.

        align mol1 & name CA, mol2 & name CA, object=aln_mol1_to_mol2
        rmsd_b mol1 & aln_mol1_to_mol2, mol2 & aln_mol1_to_mol2, byres=1

    If you specify byres=1, then it will only calculate the RMSD for alpha-carbons (i.e.
    it modifies the selection by adding "and name CA"), but will modify the B-factors for
    all atoms in each residue to be the same as the alpha-carbons.

    If you use a selection that only includes only some atoms, only those atoms will have
    their B-factors modified (even with byres=1 specified).

  """

# setting byres to true only makes sense if the selection is only for the alpha-carbons
  byres = int(byres)

# initialize the arrays used
  models = []
  coord_array = []
  chain=[]
  resi=[]
  resn=[]
  name=[]


# get models into models array
  models.append(cmd.get_model(selection1))
  models.append(cmd.get_model(selection2))

# extract coordinates and atom identification information out of the models
# loop over the states
  for i in range(len(models)):
    coord_array.append([])
    chain.append([])
    resi.append([])
    resn.append([])
    name.append([])

    if byres:
# loop over the atoms in each state
      for j in range(len(models[0].atom)):
        atom = models[i].atom[j]
        if atom.name == 'CA':
          coord_array[i].append(atom.coord)
          chain[i].append(atom.chain)
          resi[i].append(atom.resi)
          resn[i].append(atom.resn)
          name[i].append(atom.name)
    else:
# loop over the atoms in each state
      for j in range(len(models[0].atom)):
        atom = models[i].atom[j]
        coord_array[i].append(atom.coord)
        chain[i].append(atom.chain)
        resi[i].append(atom.resi)
        resn[i].append(atom.resn)
        name[i].append(atom.name)

# initialize interatomic distance array
  dist = []
  distx2 = []

# calculate the square of the fluctuation
# = the sum of the distance between the atoms in each state from the reference state, squared
  for j in range(len(coord_array[0])):
    d2 = distance2(coord_array[0][j],coord_array[1][j])
    dist.append(math.sqrt(d2))
    distx2.append(d2)

# reset all the B-factors to zero
  cmd.alter(selection1,"b=0")

  b_dict = {}
  sum_distx2 = 0

# if byres, alter all B-factors in a residue to be the same
  if byres:
    def b_lookup(chain, resi, name):
      if chain in b_dict:
        b = b_dict[chain][resi]
      else:
        b = b_dict[''][resi]
      return b
    stored.b = b_lookup

    for i in range(len(dist)):
      sum_distx2 += distx2[i]
      b_dict.setdefault(chain[0][i], {})[resi[0][i]] = dist[i]
    mean_rmsd = math.sqrt(sum_distx2/len(dist))
    cmd.alter(selection1,'%s=stored.b(chain,resi,name)' % ('b'))

# if not byres, alter all B-factors individually according to the RMSD for each atom
  else:
    quiet=0
    def b_lookup(chain, resi, name):
      if chain in b_dict:
        b = b_dict[chain][resi][name]
      else:
        b = b_dict[''][resi][name]
      return b
    stored.b = b_lookup

    for i in range(len(dist)):
      sum_distx2 += distx2[i]
      try:
        b_dict.setdefault(chain[0][i], {}).setdefault(resi[0][i], {})[name[0][i]] = dist[i]
      except KeyError:
        print(chain[0][i],resi[0][i],name[0][i])
    mean_rmsd = math.sqrt(sum_distx2/len(dist))
    cmd.alter(selection1,'%s=stored.b(chain,resi,name)' % ('b'))

  print("Mean RMSD for selection: %s = %g" % (selection1, mean_rmsd))

cmd.extend("rmsd_b",rmsd_b)
