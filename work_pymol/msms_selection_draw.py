#! /usr/bin/env python3
# Copyright (c) 2008 Robert L. Campbell (rlc1@queensu.ca)
#
import sys

def get_color(x,inside_color=(1,1,1),outside_color=(0.5,0.5,0.5)):
# x = 1 (inside selection) , x = 0 == (outside selection), 
  if x == 0:
    color_value = outside_color
  else:
    color_value = inside_color

#    if not color_value: print x
  return color_value

def msms_obj(vert_name,face_name,draw,inside_color_tuple,outside_color_tuple):
#infile = sys.stdin
  vert_name = vert_name.replace('.vert','')
  vert_infile = open(vert_name + '.vert')
# strip off header lines
  vert_line1 = vert_infile.readline()
  vert_line2 = vert_infile.readline()
  vert_line3 = vert_infile.readline()

# Vertices file contains:
# vert1,vert2,vert3,norm1,norm2,norm3,face_num,sphere_num,toric_flag,atom_ID
  vertices = []
  for line in vert_infile.readlines():
    data = line.strip().split()
    vertex = [float(x) for x in data[0:3]]
    normal = [float(x) for x in data[3:6]]
    face_number = int(data[6])
    sphere_index = int(data[7])
    toric_flag = int(data[8])
    atom_name = data[9]
    try:
      draw_flag = float(data[10])
    except IndexError:
      draw_flag = 1
    vertices.append((vertex,normal,face_number,sphere_index,toric_flag,atom_name,draw_flag))

  face_name = face_name.replace('.face','')
  face_infile = open(face_name + '.face')
# strip off header lines
  face_line1 = face_infile.readline()
  face_line2 = face_infile.readline()
  face_line3 = face_infile.readline()

# Face file contains triangles composed of vertices (1-based index):
# vert_index1,vert_index2,vert_index3,toric_flag,face_num(of analytical surf)
  faces = []
  for line in face_infile.readlines():
    data = [int(x) for x in line.strip().split()]
    faces.append(data)


  msms_obj = []
  msms_obj.extend([BEGIN,TRIANGLES])

  if draw == 1:
    for face in faces:
      if vertices[face[0]-1][6] != 0:
# add 3 vertex coordinates
# need to subtract one from face[0] to get index into vertices array
        msms_obj.extend([VERTEX])
        msms_obj.extend(vertices[face[0]-1][0])
        msms_obj.extend([NORMAL])
        msms_obj.extend(vertices[face[0]-1][1])
        msms_obj.extend([COLOR])
        msms_obj.extend(inside_color_tuple)

        msms_obj.extend([VERTEX])
        msms_obj.extend(vertices[face[1]-1][0])
        msms_obj.extend([NORMAL])
        msms_obj.extend(vertices[face[1]-1][1])
        msms_obj.extend([COLOR])
        msms_obj.extend(inside_color_tuple)

        msms_obj.extend([VERTEX])
        msms_obj.extend(vertices[face[2]-1][0])
        msms_obj.extend([NORMAL])
        msms_obj.extend(vertices[face[2]-1][1])
        msms_obj.extend([COLOR])
        msms_obj.extend(inside_color_tuple)
  elif draw == 2:
    for face in faces:
# add 3 vertex coordinates
# need to subtract one from face[0] to get index into vertices array
      msms_obj.extend([VERTEX])
      msms_obj.extend(vertices[face[0]-1][0])
      msms_obj.extend([NORMAL])
      msms_obj.extend(vertices[face[0]-1][1])
# use get_color function to get either the "inside" or "outside" colour, depending on the value of "vertices[face[i]-1][6]
      color = get_color(vertices[face[0]-1][6],inside_color_tuple,outside_color_tuple)
      msms_obj.extend([COLOR])
      msms_obj.extend(color)

      msms_obj.extend([VERTEX])
      msms_obj.extend(vertices[face[1]-1][0])
      msms_obj.extend([NORMAL])
      msms_obj.extend(vertices[face[1]-1][1])
      color = get_color(vertices[face[1]-1][6],inside_color_tuple,outside_color_tuple)
      msms_obj.extend([COLOR])
      msms_obj.extend(color)

      msms_obj.extend([VERTEX])
      msms_obj.extend(vertices[face[2]-1][0])
      msms_obj.extend([NORMAL])
      msms_obj.extend(vertices[face[2]-1][1])
      color = get_color(vertices[face[2]-1][6],inside_color_tuple,outside_color_tuple)
      msms_obj.extend([COLOR])
      msms_obj.extend(color)
  msms_obj.extend([END])
  return msms_obj

# if running as a standalone program
if __name__ == '__main__':
  vert_name = sys.argv[1]
  face_name = sys.argv[2]

# figure out a way to easily pass in these values for running
# without pymol
  BEGIN='BEGIN'
  TRIANGLES='TRIANGLES'
  VERTEX='VERTEX'
  NORMAL='NORMAL'
  COLOR='COLOR'
  END='END'

  msms_obj = msms_obj(vert_name,face_name)
  color = [float(x) for x in color]
  print(msms_obj[0],msms_obj[1])
  for i in range(2,len(msms_obj)-1,12):
    for j in range(0,11):
      print(msms_obj[i+j], end=' ')
    print(msms_obj[i+11])
  print(msms_obj[-1])

# if not calling directly, then extend the def for use in PyMOL
else:
  from pymol.cgo import *    # get constants
  from pymol import cmd

  def msms_selection_draw(vert_name,face_name,obj_name='',draw=1,inside_color='white',outside_color='gray'):
#  def msms_selection_draw(vert_name,face_name,obj_name=''):
    """
    msms_selection_draw vert_name, face_name, [obj_name ], [draw], [inside_color],[outside_color]
    """
    draw = int(draw)
    if '(' in inside_color:
      print("Using inside_color as tuple:", inside_color.replace('(','').replace(')','').split(','))
      inside_color_tuple = list(map(float, inside_color.replace('(','').replace(')','').split(',')))
    else:
      inside_color_tuple = cmd.get_color_tuple(inside_color)

    if '(' in outside_color:
      print("Using outside_color as tuple:", outside_color.replace('(','').replace(')','').split(','))
      outside_color_tuple = list(map(float, outside_color.replace('(','').replace(')','').split(',')))
    else:
      outside_color_tuple = cmd.get_color_tuple(outside_color)

#    print inside_color_tuple,outside_color_tuple
    msms_cgo = msms_obj(vert_name,face_name,draw,inside_color_tuple,outside_color_tuple)
    if obj_name == '':
      obj_name = vert_name + 'msms'
    cmd.load_cgo(msms_cgo,obj_name)

  cmd.extend('msms_selection_draw',msms_selection_draw)
