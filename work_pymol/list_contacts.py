# /usr/bin/env python3
# Copyright (c) 2010 Robert L. Campbell
from pymol import cmd

def list_contacts(selection1,selection2,cutoff=4,contact_name='contacts'):
  """
  USAGE

  list_contacts selection1, selection2, [cutoff (default=4], [contact_name (default='contacts']

  e.g.
    list_contacts 1abc & c. a &! r. hoh, 1xyz & c. b &! r. hoh, cutoff=3.2, contact_name=abc-xyz_contacts
  """
  cutoff=float(cutoff)
  contacts = cmd.find_pairs(selection1,selection2,mode=0,cutoff=cutoff)
# sort the list for easier reading
  contacts.sort(key=lambda x:x[0][1])

  for pairs in contacts:
    cmd.iterate("%s and index %s" % (pairs[0][0],pairs[0][1]), 'print("%s/%3s`%s/%s " % (chain,resn,resi,name),end=" ")')
    cmd.iterate("%s and index %s" % (pairs[1][0],pairs[1][1]), 'print("%s/%3s`%s/%s " % (chain,resn,resi,name),end=" ")')
    print("%.2f" % cmd.distance(contact_name,"%s and index %s" % (pairs[0][0],pairs[0][1]),"%s and index %s" % (pairs[1][0],pairs[1][1])))

cmd.extend("list_contacts",list_contacts)
