#! /usr/bin/env python3
# Copyright (c) 2010 Robert L. Campbell (rlc1@queensu.ca)
from pymol.cgo import *    # get constants
from pymol import cmd

def draw_sphere_cgo(name,coor,radius=4.0,color=(1,1,1)):
  radius = float(radius)
  if type(color) == type(''):
    color = list(map(float,color.replace('(','').replace(')','').split(',')))
  x = list(map(float,coor))

  cgo_sphere = [COLOR, color[0],color[1],color[2],SPHERE,x[0],x[1],x[2],radius]

  cmd.load_cgo(cgo_sphere,name)

def draw_sphere(name, selection='(pk1)',radius=4,color=(1,1,1)):
  coor = cmd.get_model(selection).atom[0].coord
  draw_sphere_cgo(name, coor,radius,color)

cmd.extend('draw_sphere',draw_sphere)
