#! /usr/bin/python3
# Copyright (c) 2008 Robert L. Campbell (rlc1@queensu.ca)
# with corrections and additions by Steve Darnell
#
import os,time,sys,re
from pymol import cmd

def convert_pdb_to_xyzrn(pdb_file,xyzrn_file,atmnumbers_file):
  amino_acid_names_list = ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU',
                           'MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']
  npats = 0
  explicit_rad = {}
  united_rad = {}
  respat = {}
  atmpat = {}
  atmnum = {}
  numfile = atmnumbers_file
  lines = open(numfile).readlines()

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

  pdb_lines = open(pdb_file).readlines()
  xyzrn = open(xyzrn_file,'w')
  for line in pdb_lines:
    if line[0:4] == 'ATOM' or line[0:4] == 'HETA' or line[0:4] == 'atom' or line[0:4] == 'heta':
      x = float(line[30:38]) # 31 - 38        Real(8.3)     x
      y = float(line[38:46]) # 39 - 46        Real(8.3)     y
      z = float(line[46:54]) # 47 - 54        Real(8.3)     z
      resname = line[17:20]  # 18 - 20        Residue name  resName
      aname = line[12:16]    # 13 - 16        Atom          name
      atype = aname

# kludge to look for [0-9](H|D) in atom name for Hydrogens
      if h_match.match(aname[0:2]):
        atype = 'H'
# kludge to look for H... in atom name for non-conforming Hydrogens of amino acids
      if aname[0:1] == 'H' and resname in amino_acid_names_list:
        atype = 'H'

      chain = line[21]          # 22             Character     chainID
      resnum = int(line[22:26]) # 23 - 26        Integer       resSeq
      resname = resname.strip()
      atype = atype.strip()

      for pat in range(npats):
        if re.search(atmpat[pat],atype) and re.search(respat[pat],resname):
#        print("Found", pat, atmpat[pat],atype,respat[pat],resname)
          break
      #print(atmpat[pat],atype)
      #print(respat[pat],resname)
      if pat == npats:
#not found
        sys.stderr.write("pdb_to_xyzr: error, file %s line %s residue %d atom pattern %s %s was not found in %s\n" %( filename,NR,resnum,resname,aname,numfile) )
        print(x,y,z,0.01)
      else:
        xyzrn.write("%f %f %f %f %d %s_%s_%d_%s\n" % (x,y,z,united_rad[atmnum[pat]],1,chain,resname,resnum,aname.strip()))
  xyzrn.close()

def calc_msms(object_name,probe=1.4,density=1,keep_tmp=0,debug=0):
  probe=float(probe)
  density=float(density)
  tmp_name = "tmp_%s_%s" % (object_name,str(time.time()).replace('.','_'))
  tmp_pdb_file = tmp_name + ".pdb"
  tmp_xyzrn_file = tmp_name + ".xyzrn"
  tmp_msms_area_file = tmp_name + ".area"
  tmp_msms_vert_file = tmp_name + ".vert"
  tmp_msms_face_file = tmp_name + ".face"
  cmd.save(tmp_pdb_file,object_name,1,'pdb')
  path_list = os.environ['PATH'].split(os.pathsep)
  path_list.append('C:\Program Files\MSMS')
# users may need to append another path specific to their system here, for example:
#  path_list.append('C:\Users\Rob\bin\msms_win32_2.6.1')
  msms_path = ''
  atmnumbers_file = ''
  found_msms = 0
  found_atmtypenumbers = 0

  for p in path_list:
   if os.path.exists("%s/msms" % p):
     found_msms = 1
     msms_path = p+"/msms"
   elif os.path.exists("%s/msms.exe" % p):
     found_msms = 1
     msms_path = p+"/msms.exe"
   if os.path.exists("%s/atmtypenumbers" % p):
     found_atmtypenumbers = 1
     atmnumbers_file = p+"/atmtypenumbers"

  if not found_msms:
    print("I could not find the msms program.")
    print("You must install and/or set up msms before starting PyMOL")
    print("Path contains:")
    for p in path_list:
      print(p)
    return None,None,None,None

  elif not found_atmtypenumbers:
    print("I could not find the atmtypenumbers file")
    print("You must install and/or set up msms before starting PyMOL")
    print("Path contains:")
    for p in path_list:
      print(p)
    return None,None,None,None

  else:
    convert_pdb_to_xyzrn(tmp_pdb_file,tmp_xyzrn_file,atmnumbers_file)
# Just in case someone has put the executable in a directory with spaces in the names
# (likely the case in Windows) enclose msms_path in escaped double quotes
    msms_path = '\"%s\"' % (msms_path,)
    msms_cmd = '%s -probe_radius %f -density %f -if %s -af %s -of %s' % (
                msms_path,probe,density,tmp_xyzrn_file,tmp_msms_area_file,tmp_msms_vert_file.replace('.vert',''))
    if debug:
      print("MSMS_CMD: ",msms_cmd)
    msms_output = os.popen(msms_cmd).readlines()
    for line in msms_output:
      data = line.split()
      if len(data) == 3 and data[1] == 'ses_volume:':
        volume = data[2]
    try:
      msms_area = open(tmp_msms_area_file).readlines()
      if not keep_tmp:
        os.unlink(tmp_pdb_file)
        os.unlink(tmp_msms_area_file)
        os.unlink(tmp_xyzrn_file)
      else:
        print("MSMS surface data is in \n%s \n%s \n%s" % (tmp_msms_area_file,tmp_msms_vert_file,tmp_msms_face_file))
        print("You could display the surface with the msms_cgo.py script")
      if debug:
        for line in msms_output:
          print(line[:-1])
    except:
      print("There seems to be a problem with the msms output\n")
      for line in msms_output:
        print(line[:-1])
      #sys.exit(1)
      raise Exception("ERROR: MSMS failed for "+object_name)

    area_dict = get_area_dict(msms_area)
    return volume,area_dict,tmp_name,tmp_msms_vert_file,tmp_msms_face_file

def get_area_dict(msms_area):
  area_dict = {}
  for line in msms_area[1:]:
    atnum,ses,sas,name = line.strip().split()
    atnum = int(atnum)
    ses = float(ses)
    sas = float(sas)
    area_dict[name] = [ses,sas]
  return area_dict

def get_selection_list(object_name,selection):

  selection = "%s & %s" % (object_name,selection)
# exclude waters
  waters = 'r. hoh+wat'
  selection = '%s &! %s' % (selection,waters)
# or could do:
#  amino_acids = 'r. ala+arg+asn+asp+cys+gln+glu+gly+his+ile+leu+lys+met+phe+pro+ser+thr+trp+tyr+val'
#  ions = 'r. ca+na+cl+k'
#  selection = '%s & %s' % (selection,amino_acids,ions)

  m = cmd.get_model(selection)
  selection_list = []
  for i in m.atom:
    selection_list.append('%s_%s_%s_%s' % (i.chain,i.resn,i.resi,i.name))

  return selection_list

def atom_type(name_string):

  amino_acid_names = ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU',
                           'MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']
  hydrophobics = {'*': ['c*']}
  negative_amino_acids = ['ASP','GLU']
  positive_amino_acids = ['ARG','LYS','HIS']

#  charged_negative = {'GLU':['OE1','OE2'],'ASP':['OD1','OD2']}
#  charged_positive = {'ARG':['NH1','NH2','NE'],'LYS':['NZ']}

  values = name_string.split('_')
  if len(values) == 4:
    chain,resn,resi,name = values
  else:
    resn,resi,name = values
    chain = ''

  if resn in amino_acid_names and name[0] == 'C':
    atom_type = 'hydrophobic'
  elif resn in positive_amino_acids and name in ['NE','NH1','NH2','NZ']:
    atom_type = 'positive'
  elif resn in negative_amino_acids and name in ['OD1','OD2','OE1','OE2']:
    atom_type = 'negative'
  elif resn in amino_acid_names and name[0] == 'O' :
    atom_type = 'polar_O'
  elif resn in amino_acid_names and name[0] == 'N':
    atom_type = 'polar_N'
  elif resn in ['CYS','MET'] and name[0] == 'S':
    atom_type = 'sulfur'
  else:
    atom_type = 'other'
  return atom_type

def calc_msms_area(object_name,selection='',file_name='',probe=1.4,density=1.0,calc_residue_area=0,keep=0,debug=0,draw='',draw_object='',inside_color='white',outside_color='gray'):
  """

  AUTHOR

    Robert L. Campbell

  USAGE

    calc_msms_area object_name='object', selection='',file_name='',probe=1.4,density=1.0,
               calc_residue_area=0, keep=0, debug=0, draw='',
               inside_color='white', outside_color='gray'

    MSMS will be used to calculate the area and volume for the 'object',
    and this script calculates the total SES (solvent-excluded surface) and SAS
    (solvent-accessible surface) areas for just the selection.  If the
    'selection' is omitted, then the whole 'object' will be used.

    e.g.  calc_msms_area 1dvi, 1dvi & c. a within 5 of 1dvi & c. b, probe=1.4, draw=part

    It will also report the total SES volume of the object, regardless of
    what selection is specified.

    If you input a file_name, it will read the areas from that, rather than
    recalculating them.  You can change the probe radius (default=1.5A) and
    probe density (default = 1) values that MSMS uses, by setting them on
    the command line.

    If the "calc_residue_area" option is set to 1, then a file will be
    created that summarizes the surface area per residue.

    The "draw" option allows the generation of a CGO object showing
    the surface calculated.  Setting draw to 1 or part shows only the
    area calculated (coloured by "inside_color"), while setting draw to
    2 or full will draw the complete surface for the object with the
    area of interest colored by "inside_color" and the rest coloured by
    "outside_color".

    The variable "draw_object" can be used to specify the name of the
    msms cgo object drawn if draw > 0.

    The "keep" option does not delete the temporary xyzrn and area files
    that MSMS generates.

    The "debug" option prints the MSMS log file to the terminal.

  """
# put in option to feed msms area file directly, rather than calculating on the fly (for windows?)
# also look at possibility of bypassing pdb_to_xyzrn, but I don't think so
  probe=float(probe)
  density=float(density)
  keep=int(keep)
  debug=int(debug)
  calc_residue_area = int(calc_residue_area)
  draw_cgo = 0
  if draw != '' or draw == '0':
    if draw == 'part' or draw == '1':
      draw_cgo = 1
    elif draw == 'full' or draw == '2':
      draw_cgo = 2
  if file_name != '':
    area_dict = get_area_dict(open(file_name).readlines())
    volume = 'Unknown'
  else:
    volume,area_dict,tmp_name,tmp_msms_vert_file,tmp_msms_face_file = calc_msms(object_name,probe,density,keep,debug)
  if area_dict:
    if selection == '':
      selection = object_name

    selection_list = get_selection_list(object_name,selection)
    area_ses = 0
    area_sas = 0
    if debug:
      print("Length of selection_list", len(selection_list))

    sumarea = {'hydrophobic':0,'positive':0,'negative':0,'polar_O':0,'polar_N':0,'sulfur':0,'other':0}
    for name in selection_list:
      if debug: print(name,area_dict[name][0],area_dict[name][1])
      if name in list(area_dict.keys()):
        area_ses = area_ses + area_dict[name][0]
        area_sas = area_sas + area_dict[name][1]
        area_type = atom_type(name)
        sumarea[area_type] += area_dict[name][0]

    if draw_cgo > 0:
      new_msms_vert_file = tmp_msms_vert_file.replace('tmp','new')
      new_vert = open(new_msms_vert_file,'w')
      lines = open(tmp_msms_vert_file).readlines()
      for line in lines[0:3]:
        new_vert.write("%s" % line)
      for line in lines[3:]:
        if line.split()[9] in selection_list:
          new_vert.write("%s %d\n" % (line[:-1],1))
        else:
          new_vert.write("%s %d\n" % (line[:-1],0))
      new_vert.close()
      if not draw_object:
        obj_name = object_name + '_msms'
      else:
        obj_name = draw_object
# import this here, as it is unnecessary otherwise.
      import msms_selection_draw

      msms_selection_draw.msms_selection_draw(new_msms_vert_file,tmp_msms_face_file,obj_name,draw=draw_cgo,inside_color=inside_color,outside_color=outside_color)
      print("Creating cgo file for selection: %s " % (selection))
      if not keep:
#        print("closing temporary file:",new_msms_vert_file)
        os.unlink(new_msms_vert_file)

    if not keep:
#      print("closing temporary files:",tmp_msms_vert_file,tmp_msms_face_file)
      os.unlink(tmp_msms_vert_file)
      os.unlink(tmp_msms_face_file)

    if calc_residue_area:
      import msms_residue_area
      msms_residue_area.calc_residue_area(area_dict,filename=tmp_name+"_residue.area")

    print("For the selection: %s " % (selection))
    print("Total SES (solvent excluded surface) area is %8.2f" % (area_ses))
    print("Total SAS (solvent accessible surface) area is %8.2f" % (area_sas))
    print("\nFor the object: %s " % (object_name))
    print("Total SES volume is %s" % (volume))
    print("\nSES Areas by atom-type and percentage of total")
    for at in list(sumarea.keys()):
      print("%11s: %10.3f %5.1f%%" % (at,sumarea[at],sumarea[at]*100/area_ses))
    print("")
    return (area_ses, area_sas)

def calc_msms_area_for_directory(logfile, dir='.', ext='.pdb'):
  '''
  Calculate the MSMS surface for each PDB file in the given directory and
  records the molecular surface area and solvent accessible surface area for
  each file.

  logfile  Record surface areas to file
  dir      Directory containing PDB files
  ext      File extension for PDB files
  '''
  log = open(logfile,'w')
  files = os.listdir(dir)
  for file in files:
    if file.endswith(ext):
      cmd.load(os.path.join(dir,file))
      objname = file[:file.find(ext)]
      try:
        molec_surf_area, solvent_acc_surf_area = calc_msms_area(objname)
      except Exception as e:
        print(e)
        molec_surf_area, solvent_acc_surf_area = (float('nan'),float('nan'))
      cmd.delete(objname)
      log.write('%s\t%.2f\t%.2f\n' % (objname, molec_surf_area, solvent_acc_surf_area))
  log.close()

cmd.extend("calc_msms_area",calc_msms_area)
cmd.extend("calc_msms_area_for_directory", calc_msms_area_for_directory)
