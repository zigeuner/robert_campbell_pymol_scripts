#! /usr/bin/env python
# Copyright (c) 2005 Robert L. Campbell

from pymol import cmd, stored

def fitting_v1(obj1,select1,obj2,select2):
  """
DESCRIPTION

  "fitting_v1" allows the superpositioning of object1 onto object2 using
  the atoms in selection1 and selection2 (side chains are ignored).
  For both selections, the residue names are changed to 'GLY',residue
  numbers are set to be sequential starting at '1', chain identifiers
  are set to '', segment identifiers are set to '1fit', and alt ids are
  set to '', temporarily.  This allows the normal "fit" command to work.
  They are reset after "fit" is run and two new objects are created
  showing the selected and fit atoms.

  It is important that the beginning residue numbers specify the
  aligned residues, but the ending numbers are not critical.  The
  shorter of the two selections is used in the fit calculation.

USAGE 

  fitting_v1 object1, selection1, object2, selection2


EXAMPLES

  fitting_v1 1xuu, c. a & i. 296-309, 1ame, i. 8-21

  """

  list_m = []
  list_n = []

  backbone = 'n. n+ca+c+o &! r. hoh+wat'
  select1 += ' & %s' % backbone
  select2 += ' & %s' % backbone
  m=cmd.get_model("%s & %s" % (obj1,select1))
  n=cmd.get_model("%s & %s" % (obj2,select2))

# for the atoms to be used in fit:
# store id, chain, resn, resi, name, segi, alt
  for at in m.atom:
    list_m.append((at.id,at.chain,at.resn,at.resi,at.name,at.segi, at.alt))
  for at in n.atom:
    list_n.append((at.id,at.chain,at.resn,at.resi,at.name,at.segi, at.alt))

  print obj1, "and", select1,"to", obj2, "and", select2

# set a new segi for the atoms to be used in fit command and to allow resetting later
  seg_fit="1fit"

# set a counter for residue numbers and increment resi_count when the residue number is different from the previous one.
  last_chain = ''
  last_resn = ''
  last_resi = -9999
  resi_count = 0
# Alter chain,resn,resi & segi
  for atoms_m in list_m:
    if atoms_m[1] != last_chain or atoms_m[2] != last_resn or atoms_m[3] != last_resi:
      last_chain = atoms_m[1]
      last_resn = atoms_m[2]
      last_resi = atoms_m[3]
      resi_count += 1
    cmd.do("alter %s & id %s, chain=''" % (obj1,atoms_m[0]))
    cmd.do("alter %s & id %s, resn='GLY'" % (obj1,atoms_m[0]))
    cmd.do("alter %s & id %s, resi=%s" % (obj1,atoms_m[0],str(resi_count)))
    cmd.do("alter %s & id %s, segi='%s'" % (obj1,atoms_m[0],seg_fit))
    cmd.do("alter %s & id %s, alt=''" % (obj1,atoms_m[0]))
  num_m = resi_count

  last_chain = ''
  last_resn = ''
  last_resi = -9999
  resi_count = 0
  for atoms_n in list_n:
    if atoms_n[1] != last_chain or atoms_n[2] != last_resn or atoms_n[3] != last_resi:
      last_chain = atoms_n[1]
      last_resn = atoms_n[2]
      last_resi = atoms_n[3]
      resi_count += 1
    cmd.do("alter %s & id %s, chain=''" % (obj2,atoms_n[0]))
    cmd.do("alter %s & id %s, resn='GLY'" % (obj2,atoms_n[0]))
    cmd.do("alter %s & id %s, resi=%s" % (obj2,atoms_n[0],str(resi_count)))
    cmd.do("alter %s & id %s, segi='%s'" % (obj2,atoms_n[0],seg_fit))
    cmd.do("alter %s & id %s, alt=''" % (obj2,atoms_n[0]))
  num_n = resi_count

  print "%s & resi %s-%s & segi %s & %s\n" % (obj1,'0',str(num_m),seg_fit,backbone),
  print "%s & resi %s-%s & segi %s & %s\n" % (obj2,'0',str(num_n),seg_fit,backbone) 
  rms = cmd.fit("%s & resi %s-%s & segi %s & %s" % (obj1,'0',str(num_m),seg_fit,backbone),
                "%s & resi %s-%s & segi %s & %s" % (obj2,'0',str(num_n),seg_fit,backbone) ,quiet=0)

  cmd.delete("%s_fitv1" % obj1)
  cmd.delete("%s_fitv1" % obj2)
# create new objects to show the fit atoms
  cmd.create("%s_fitv1" % obj1, "%s & resi %s-%s & segi %s & %s" % (obj1,'0',str(num_m),seg_fit,backbone))
  cmd.create("%s_fitv1" % obj2, "%s & resi %s-%s & segi %s & %s" % (obj2,'0',str(num_n),seg_fit,backbone))

# reset chain,resn,resi & segi from stored lists
  for atoms_m in list_m:
    cmd.do("alter %s & id %s, chain='%s'" % (obj1,atoms_m[0],atoms_m[1]))
    cmd.do("alter %s & id %s, resn='%s'" % (obj1,atoms_m[0],atoms_m[2]))
    cmd.do("alter %s & id %s, resi=%s" % (obj1,atoms_m[0],atoms_m[3]))
    cmd.do("alter %s & id %s, segi='%s'" % (obj1,atoms_m[0],atoms_m[5]))
    cmd.do("alter %s & id %s, alt='%s'" % (obj1,atoms_m[0],atoms_m[6]))
  for atoms_n in list_n:
    print atoms_n
    cmd.do("alter %s & id %s, chain='%s'" % (obj2,atoms_n[0],atoms_n[1]))
    cmd.do("alter %s & id %s, resn='%s'" % (obj2,atoms_n[0],atoms_n[2]))
    cmd.do("alter %s & id %s, resi=%s" % (obj2,atoms_n[0],atoms_n[3]))
    cmd.do("alter %s & id %s, segi='%s'" % (obj2,atoms_n[0],atoms_n[5]))
    cmd.do("alter %s & id %s, alt='%s'" % (obj2,atoms_n[0],atoms_n[6]))

  print "RMSD for fitting selection %s of %s onto \n                 selection %s of %s = %6.3f" % (select1, obj1, select2, obj2, rms)

cmd.extend("fitting_v1",fitting_v1)
