#! /usr/bin/python
# Copyright (c) 2004 Robert L. Campbell

from pymol.cgo import *
from pymol import cmd

from numpy import *

# get cell, wavelength info from a separate file
#from cell_info import *
##########################################################################
from cctbx import uctbx, sgtbx

a = 30.
b = 30.
c = 40.
alpha = 90.
beta = 90.
gamma = 90.
sg = 'P422'

dmin = 2.
# set angular position of starting point of ewald sphere
phi0 = 0.
phi0_rad = phi0 * pi / 180.

# expand by 10 ease of viewing in PyMOL
dmin_recip = 10. / dmin
dminsq_recip = dmin_recip * dmin_recip

wavelength0 = .2
wavelength1 = 2.5
r0 = 10./wavelength0
r1 = 10./wavelength1
r0sq = r0*r0
r1sq = r1*r1

SgSymb = sgtbx.SpaceGroupSymbols(sg)
SgOps = sgtbx.SpaceGroup(SgSymb)
#print SgOps

ASU = sgtbx.ReciprocalSpaceASU(SgOps.Info())

U = uctbx.UnitCell((a,b,c,alpha,beta,gamma))
U_recip = uctbx.UnitCell(U.getParameters(1))

(hmax,kmax,lmax) = U.MaxMillerIndices(dmin)
print "(Hmax, Kmax, Lmax) ", (hmax,kmax,lmax)

phirotmat0 = array([[cos(phi0_rad), 0, sin(phi0_rad)],
                    [ 0, 1, 0],
                    [-sin(phi0_rad), 0, cos(phi0_rad)]])
#phirotmatsweep = array([[cos(phisweep_rad), 0, sin(phisweep_rad)],
#                        [ 0, 1, 0],
#            [-sin(phisweep_rad), 0, cos(phisweep_rad)]])
#
c0 = dot(phirotmat0,array([0,0,-r0]))
c1 = dot(phirotmat0,array([0,0,-r1]))
#c1 = dot(phirotmatsweep,c0)

#c0 = array([0,0,-r0])
#c1 = array([0,0,-r1])

print "c0 = ", c0, "    c1 = ", c1

# check whether second ewald sphere is beyond 180-cusp angle away from 1st
#if phisweep > 180 - 2 * asin(wavelength/(2*dmin)


def distsq(a,b):
  distsq = sum((a-b)*(a-b))
  #dist = sqrt(sum((a-b)*(a-b)))
  return distsq

def dist(a,b):
  dist = sqrt(sum((a-b)*(a-b)))
  return dist

##########################################################################

def make_cross(x,y,z):
  dx = r1/100.
  str = [VERTEX, x-dx, y, z, VERTEX, x+dx, y, z,
      VERTEX, x, y-dx, z, VERTEX, x, y+dx, z,
      VERTEX, x, y, z-dx, VERTEX, x, y, z+dx, ]
  return str

def within_sph1(x,y,z,xsq,ysq,zsq,dsq):
  try:
    if distsq(array([x,y,z]),c0) < r0sq and distsq(array([x,y,z]),c1) > r1sq:
      return 1
    else:
      return 0
  except ValueError, err:
    print "x, y, z, d, r,  xsq, ysq, zsq, dsq, r0sq "
    print x,y,z,1/sqrt(dsq),1/r,xsq,ysq,zsq,dsq,r0sq
    print "\n\n"
    print err
    return 0

def within_sph2(x,y,z,xsq,ysq,zsq,dsq):
  try:
    if distsq(array([x,y,z]),c0) < rsq and distsq(array([x,y,z]),c1) > rsq:
      return 1
    else:
      return 0
  except ValueError, err:
    print "x, y, z, d, r,  xsq, ysq, zsq, dsq, rsq "
    print x,y,z,1/sqrt(dsq),1/r,xsq,ysq,zsq,dsq,rsq
    print "\n\n"
    print err
    return 0

def draw_ewald0():
  ewald0_obj = [COLOR, .8, .8, .8, SPHERE, c0[0],c0[1],c0[2],r0]
  cmd.load_cgo(ewald0_obj,'ewald0')

def draw_ewald1():
  ewald1_obj = [COLOR, 1.0, .7, .5, SPHERE, c1[0],c1[1],c1[2],r1 ]
  cmd.load_cgo(ewald1_obj,'ewald1')

def draw_limsph():
  limsph_obj = [COLOR, .0, .5, .5, SPHERE, 0,0,0,dmin_recip]
  cmd.load_cgo(limsph_obj,'limsph')


def draw_axes():

  axes_obj = [BEGIN, LINES, 
        COLOR, 1., 1., 1.,
        VERTEX, 0., 0., 0.,
        COLOR, 1., 0., 0.,
        VERTEX, r, 0., 0.,

        COLOR, 1., 1., 1.,
        VERTEX, 0., 0., 0.,
        COLOR, 0., 1., 0.,
        VERTEX, 0., r, 0.,

        COLOR, 1., 1., 1.,
        VERTEX, 0., 0., 0.,
        COLOR, 0., 0., 1.,
        VERTEX, 0., 0., r,
        END,]
  cmd.load_cgo(axes_obj,'axes')

def draw_axes_cyl():

  axes_cyl_obj = [
       CYLINDER, 0., 0., 0., r, 0., 0., r1/100, 1.0, 0.7, 0.7, 1.0, 0.0, 0.,
       CYLINDER, 0., 0., 0., 0., r, 0., r1/100, 0.7, 1.0, 0.7, 0., 1.0, 0.,
       CYLINDER, 0., 0., 0., 0., 0., r, r1/100, 0.7, 0.7, 1.0, 0., 0.0, 1.0,
  ]
  cmd.load_cgo(axes_cyl_obj,'axes')

def draw_refl():
  
  out_obj = [BEGIN, POINTS, COLOR, 1.0, 1.0, 1.0,]

  # objects within limiting resolution sphere:
  # not collected
  in0_obj = [BEGIN, POINTS, COLOR, 1.0, 1.0, 0.0, ]
  # pass through leading edge of sphere
  in1_obj = [BEGIN, LINES, COLOR, 1.0, 0.0, 0.0, ]
  #in1_obj = [BEGIN, POINTS, COLOR, 1.0, 0.0, 0.0, ]
  # pass through trailing edge of sphere
  #in2_obj = [BEGIN, LINES, COLOR, 0.0, 1.0, 0.0, ]
  asymm_obj = [BEGIN, LINES, COLOR, 1.0, 0.5, 0.5, ]

  for h in range(-hmax,hmax):
    for k in range(-kmax,kmax):
      for l in range(-lmax,lmax):
        (x,y,z) = U_recip.orthogonalize((h,k,l))
        # expand by 10 ease of viewing in PyMOL
        (x,y,z) = (x*10.,y*10.,z*10.)
        xsq = x*x
        ysq = y*y
        zsq = z*z

        string = [VERTEX, x,y,z]

        dsq = xsq + ysq + zsq
        
        if dsq <= dminsq_recip:
          if ASU.isInASU((h,k,l)):
            asymm_obj.extend(make_cross(x,y,z))
#            asymm_obj.extend(string)
#          else:
#            in1_obj.extend(string)

          if within_sph1(x,y,z,xsq,ysq,zsq,dsq):
            in1_obj.extend(make_cross(x,y,z))
#            in1_obj.extend(string)
          else:
            in0_obj.extend(string)
        else:
          out_obj.extend(string)


  out_obj.extend([END,])
  in0_obj.extend([END,])
  in1_obj.extend([END,])
  #in2_obj.extend([END,])
  asymm_obj.extend([END,])
  
  cmd.load_cgo(out_obj,'out')
  cmd.load_cgo(in0_obj,'in0')
  cmd.load_cgo(in1_obj,'in1')
  #cmd.load_cgo(in2_obj,'in2')
  cmd.load_cgo(asymm_obj,'asymm')

#draw_axes()
draw_axes_cyl()
draw_ewald0()
draw_ewald1()
draw_limsph()
draw_refl()

