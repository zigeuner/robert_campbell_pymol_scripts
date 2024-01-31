#! /usr/bin/env python3
# Copyright (c) 2008 Robert L. Campbell (rlc1@queensu.ca)
#

import sys


def sort_by_resname(x,y):
  result = cmp(x.split('_')[0],y.split('_')[0])
  if result != 0:
    return result
  else:
    result = cmp(int(x.split('_')[2]),int(y.split('_')[2]))
    return result


def calc_residue_area(area_dict,filename=''):
  """
  area_dict is a dictionary whose keys are "chain_resn_resi_atname"
  """
  ses_res_area = {}
  sas_res_area = {}
  count = 0
  for key in list(area_dict.keys()):
    count += 1
    ses = area_dict[key][0]
    sas = area_dict[key][1]

    try:
      chain,resn,resi,atname = key.split('_')
      id = '_'.join((chain,resn,resi))
      if id in ses_res_area:
        ses_res_area[id] += ses
        sas_res_area[id] += sas
      else:
        ses_res_area[id] = ses
        sas_res_area[id] = sas

    except ValueError:
      print("Error with residue name for atom number: ",count," with id: ",key)

  resnames = list(ses_res_area.keys())
  resnames.sort(sort_by_resname)
  if filename:
    file_out = open(filename,'w')
    for id in resnames:
      chain,resnam,resnum = id.split('_')
      file_out.write('    %s %4s %4s %7.2f %7.2f\n' % (chain,resnum,resnam,ses_res_area[id],sas_res_area[id]))
  else:
    for id in resnames:
      chain,resnam,resnum = id.split('_')
      print('%s %4s %4s %7.2f %7.2f' % (chain,resnum,resnam,ses_res_area[id],sas_res_area[id]))

if __name__ == "__main__":
  if len(sys.argv) > 1:
    input = open(sys.argv[1])
  else:
    input = sys.stdin
  lines = input.readlines()

  area_dict = {}
  for line in lines[1:]:
    atnum,ses,sas,name_string = line.strip().split()
    atnum = int(atnum)
    ses = float(ses)
    sas = float(sas)
    area_dict[name_string] = [ses,sas]

  calc_residue_area(area_dict)
