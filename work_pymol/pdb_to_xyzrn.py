#! /usr/bin/env python3

import sys,re,os

amino_acid_names_list = ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU',
                         'MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']


npats = 0
explicit_rad = {}
united_rad = {}
respat = {}
atmpat = {}
atmnum = {}
if 'MSMSDIR' in os.environ:
  numfile = os.environ['MSMSDIR'] + "/atmtypenumbers"
  lines = open(numfile).readlines()
else:
  print("Need to setup access to MSMS first")
  sys.exit(1)

h_match = re.compile('[ 0-9][HhDd]')

counter = 0
for line in lines:
  counter += 1
  if line.strip() != '' and line[0] != '#':
    values = line.split()
    if values[0] == 'radius':
      n = int(values[1])
      explicit_rad[n] = float(values[3])
      if len(values) <= 4 or len(values) >= 5 and values[4][0] == '#':
        united_rad[n] = explicit_rad[n]
      else:
        united_rad[n] = float(values[4])

    else:
      respat[npats] = values[0]
      if respat[npats] == "*":
        respat[npats] = ".*"
      respat[npats] = "^%s$" % respat[npats]
      atmpat[npats] = "^%s$" % values[1]
      atmpat[npats] = atmpat[npats].replace('_',' ')
      atmnum[npats] = int(values[2])
      if atmnum[npats] not in explicit_rad:
        # the key has no radius --- complain and fake one
        sys.stderr.write("pdb_to_xyzr: error in library file %s entry %s %s %s has no corresponding radius value\n"  % (numfile, values[0],values[1],values[2]))
        explicit_rad[atmnum[npats]] = 0.01
        united_rad[atmnum[npats]] = 0.01
        
      npats += 1


#pdb_lines = open('test.pdb').readlines()
pdb_lines = sys.stdin.readlines()
for line in pdb_lines:
# included check for multiple conformation and use only the 'A' atoms if present
# the output will have the ' A' appended to the relevant atom names
  if (line[0:4] == 'ATOM' or line[0:4] == 'HETA' or line[0:4] == 'atom' or line[0:4] == 'heta') and (line[16] == ' ' or line[16] == 'A'):
    x = float(line[30:39])
    y = float(line[38:47])
    z = float(line[46:55])
    resname = line[17:21]
    aname = line[12:17]
    atype = aname

# kludge to look for [0-9](H|D) in atom name for Hydrogens
    if h_match.match(aname[0:2]):
      atype = 'H'
# kludge to look for H... in atom name for non-conforming Hydrogens of amino acids
    if aname[0:1] == 'H' and resname in amino_acid_names_list:
      atype = 'H'

    chain = line[20:22]
    resnum = int(line[22:27])
    resname = resname.strip()
    atype = atype.strip()

    for pat in range(npats):
      if re.search(atmpat[pat],atype) and re.search(respat[pat],resname):
#        print "Found", pat, atmpat[pat],atype,respat[pat],resname
        break
    #print atmpat[pat],atype
    #print respat[pat],resname
    if pat == npats:
#not found
      sys.stderr.write("pdb_to_xyzr: error, file %s line %s residue %d atom pattern %s %s was not found in %s\n" %( filename,NR,resnum,resname,aname,numfile) )
      print(x,y,z,0.01)
    else:
      print("%f %f %f %f %d %s_%s_%d_%s" % (x,y,z,united_rad[atmnum[pat]],1,chain,resname,resnum,aname.strip()))

