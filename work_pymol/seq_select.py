#! /usr/bin/python
# Copyright (c) 2004 Robert L. Campbell (rlc1@queensu.ca)
# $Id: seq_select.py,v 1.5 2004/05/14 16:37:47 rlc Exp rlc $

import re
from pymol import cmd
"""
  seq_select: create a selection in an object by searching for a particular sequence

"""
try:
# try to use seq_convert if available 
# get from http://pldserver1.biochem.queensu.ca/~rlc/work/scripts/
  from seq_convert import seq3_to_seq1

except:

  def seq3_to_seq1(seq3):
    """
    Convert array of 3-letter code sequence to 1-letter code
    Return seq1 as string
    """
    seq1 = ''
    res3 = ['---','ala','asn','asp','arg','cys','gln','glu',
             'gly','his','ile','leu','lys','met','pro','phe','ser',
             'thr','trp','tyr','val','unk',
             'ALA','ASN','ASP','ARG','CYS','GLN','GLU',
             'GLY','HIS','ILE','LEU','LYS','MET','PRO','PHE','SER',
             'THR','TRP','TYR','VAL','UNK',]
# wanted to include '^M' and '^L' in the res1 string to match 'XXX','YYY' below
#           'THR','TRP','TYR','VAL','UNK','XXX','YYY']
    res1 = '-andrcqeghilkmpfstwyvxANDRCQEGHILKMPFSTWYVX'

#force seq3 to be a list, in case it is just a single amino acid in 3-letter code
    if type(seq3) != type(list()):
      seq3 = [seq3]
    for a3 in seq3:
      #      a3 = to_upper(a3)
      a3 = a3.upper()
      # strip trailing spaces in case the three letter code came from a MyPDB Protein sequence listing
      try:
        a1 = res1[res3.index(a3.strip())]
      except ValueError as err:
        a1 = 'X'
      seq1 += a1

    return seq1


#def to_upper(a):
#  lower='abcdefghijklmnopqrstuvwxyz'
#  upper='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
#
#  b = ''
#  for i in range(len(a)):
#    try:
#      b += upper[lower.index(a[i])]
#    except:
#      b += a[i]
#  return b

def get_seq(selection):

  selection += ' & n. ca & e. c'
  m = cmd.get_model(selection)
  seq3 = []
  chains = []
  resids = []
  for i in m.atom:
    seq3.append(i.resn)
    chains.append(i.chain)
    resids.append(i.resi)

  seq1 = seq3_to_seq1(seq3)
  return seq1,chains,resids

def seq_select(obj,query,name=''):
  """
AUTHOR

  Robert L. Campbell (rlc1@post.queensu.ca)

DESCRIPTION

  "seq_select" creates a named selection of a particular sequence in 
  an object.  If the sequence occurs more than once, each instance 
  will be included in the object.

USAGE 

  seq_select object, sequence [, name ]

     if 'name' is omitted, then the selection will be named by
     concatenating the object name with 'seq'

EXAMPLES

     seq_select protein, adfg, actsite

     will create a selection of the sequence 'Ala Asp Phe Gly' in 
     the object 'protein' and will name it 'actsite'

  """
  obj_seq,obj_chain,obj_resnum = get_seq(obj)
  query_re = re.compile(query.upper())

  if name == '':
    name = obj + 'seq'

  start = []
  end = []

#iterate over matches using regexps
  matches = query_re.finditer(obj_seq)
  for match in matches:
    start.append(match.start())
    end.append(match.end())

  selection_string = '(' + obj + ' '

  if len(start) > 0:
    for j in range(len(start)):
      if j > 0:
        selection_string += ' | (' + obj + ' & ('
      else:
        selection_string += ' & ( '
      if obj_chain[start[j]] != '':
        selection_string += 'c. ' + obj_chain[start[j]] + ' & i. '
      else:
        selection_string += 'i. '
      for i in range(start[j],end[j]-1):
        selection_string += str(obj_resnum[i]) + '+'
      selection_string += str(obj_resnum[end[j]-1]) + ' ))'

    cmd.select(name,selection_string)
    print('Selecting: ', selection_string, " as ", name)

  else:
    print("Sequence %s not found in object %s" % (query, obj))

################## end of functions   ################
cmd.extend('seq_select',seq_select)
