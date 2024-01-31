#! /usr/bin/env python
# Copyright (c) 2004 Robert L. Campbell

# if not using PyMOL, then import dummy cgo and cmd modules to allow it to write information
# to stdout
try:
  from pymol.cgo import *
except:
  import cgo_testing
try:
  from pymol import cmd
except:
  import cmd_testing

from cctbx import uctbx, sgtbx

# numpy doesn't seem to work with cctbx
import numpy as N

##########################################################################

class Cell_info:

  def __init__(self,cell_param=(30.,60.,100.,90.,99.,90.),sg='p2',dmin=2.0,wavelength=1.5418,phi0=0,phisweep=-90):

    self.cell_param = cell_param
    self.sg = sg
    self.dmin = dmin
    self.wavelength = wavelength

    self.r = 10./self.wavelength
    self.rsq = self.r*self.r

#    phi0 = 0.
#    phisweep = -90.
    phi0_rad = phi0 * N.pi / 180.
    phisweep_rad = phisweep * N.pi / 180.

# expand by 10 for ease of viewing in PyMOL
    self.dmin_recip = 10. / dmin
    self.dminsq_recip = self.dmin_recip * self.dmin_recip

#    self.sg='p2'
    sg_info = sgtbx.space_group_info(self.sg)
#SgOps = sgtbx.SpaceGroup(SgSymb)
    sg_info.show_summary()
#print "Space group: ", SgOps.Info().BuildLookupSymbol()
#print SgOps

    self.asu = sg_info.reciprocal_space_asu()

    self.u_cell = uctbx.unit_cell(self.cell_param)
    self.u_cell_recip = uctbx.unit_cell((self.u_cell.reciprocal_parameters()))

    (self.hmax,self.kmax,self.lmax) = self.u_cell.max_miller_indices(dmin)
    print "(Hmax, Kmax, Lmax) ", (self.hmax,self.kmax,self.lmax)

    phirotmat0 = N.array([[N.cos(phi0_rad), 0, N.sin(phi0_rad)],
                        [ 0, 1, 0],
              [-N.sin(phi0_rad), 0, N.cos(phi0_rad)]])
    phirotmatsweep = N.array([[N.cos(phisweep_rad), 0, N.sin(phisweep_rad)],
                            [ 0, 1, 0],
                [-N.sin(phisweep_rad), 0, N.cos(phisweep_rad)]])

    self.c0 = N.dot(phirotmat0,N.array([0,0,-self.r]))

    self.c1 = N.dot(phirotmatsweep,self.c0)

    print "center of ewald sphere at position 0 = ", self.c0
    print "center of ewald sphere at position 1 = ", self.c1


def distsq(a,b):
  distsq = sum((a-b)*(a-b))
  #dist = N.sqrt(sum((a-b)*(a-b)))
  return distsq

##########################################################################

def make_cross(r,x,y,z):
  dx = r/100.
  str = [cgo.VERTEX, x-dx, y, z, cgo.VERTEX, x+dx, y, z,
      cgo.VERTEX, x, y-dx, z, cgo.VERTEX, x, y+dx, z,
      cgo.VERTEX, x, y, z-dx, cgo.VERTEX, x, y, z+dx, ]
  return str

def within_sph1(cell,x,y,z,xsq,ysq,zsq,dsq):
  try:
    if distsq(N.array([x,y,z]),cell.c0) > cell.rsq and distsq(N.array([x,y,z]),cell.c1) < cell.rsq:
      return 1
    else:
      return 0
  except ValueError, err:
    print "x, y, z, d, r,  xsq, ysq, zsq, dsq, rsq "
    print x,y,z,1/N.sqrt(dsq),1/r,xsq,ysq,zsq,dsq,cell.rsq
    print "\n\n"
    print err
    return 0

def within_sph2(cell,x,y,z,xsq,ysq,zsq,dsq):
  try:
    if distsq(N.array([x,y,z]),cell.c0) < cell.rsq and distsq(N.array([x,y,z]),cell.c1) > cell.rsq:
      return 1
    else:
      return 0
  except ValueError, err:
    print "x, y, z, d, r,  xsq, ysq, zsq, dsq, rsq "
    print x,y,z,1/N.sqrt(dsq),1/r,xsq,ysq,zsq,dsq,cell.rsq
    print "\n\n"
    print err
    return 0

def draw_ewald0(cell):
  ewald0_obj = [cgo.COLOR, .8, .8, .8, cgo.SPHERE, cell.c0[0],cell.c0[1],cell.c0[2],cell.r]
  cmd.load_cgo(ewald0_obj,'ewald0')

def draw_ewald1(cell):
  ewald1_obj = [cgo.COLOR, 1.0, .7, .5, cgo.SPHERE, cell.c1[0],cell.c1[1],cell.c1[2],cell.r ]
  cmd.load_cgo(ewald1_obj,'ewald1')

def draw_limsph(cell):
  limsph_obj = [cgo.COLOR, .0, .5, .5, cgo.SPHERE, 0,0,0,cell.dmin_recip]
  cmd.load_cgo(limsph_obj,'limsph')


def draw_axes(cell):

  axes_obj = [cgo.BEGIN, cgo.LINES,
        cgo.COLOR, 1., 1., 1.,
        cgo.VERTEX, 0., 0., 0.,
        cgo.COLOR, 1., 0., 0.,
        cgo.VERTEX, cell.r, 0., 0.,

        cgo.COLOR, 1., 1., 1.,
        cgo.VERTEX, 0., 0., 0.,
        cgo.COLOR, 0., 1., 0.,
        cgo.VERTEX, 0., cell.r, 0.,

        cgo.COLOR, 1., 1., 1.,
        cgo.VERTEX, 0., 0., 0.,
        cgo.COLOR, 0., 0., 1.,
        cgo.VERTEX, 0., 0., cell.r,
        cgo.END,]
  cmd.load_cgo(axes_obj,'axes')

def draw_axes_cyl(cell):

  axes_cyl_obj = [
       cgo.CYLINDER, 0., 0., 0., cell.r, 0., 0., cell.r/100, 1.0, 0.7, 0.7, 1.0, 0.0, 0.,
       cgo.CYLINDER, 0., 0., 0., 0., cell.r, 0., cell.r/100, 0.7, 1.0, 0.7, 0., 1.0, 0.,
       cgo.CYLINDER, 0., 0., 0., 0., 0., cell.r, cell.r/100, 0.7, 0.7, 1.0, 0., 0.0, 1.0,
  ]
  cmd.load_cgo(axes_cyl_obj,'axes')

def draw_refl(cell):
  
  out_obj = [cgo.BEGIN, cgo.POINTS, cgo.COLOR, 1.0, 1.0, 1.0,]

  # objects within limiting resolution sphere:
  # not collected
  in0_obj = [cgo.BEGIN, cgo.POINTS, cgo.COLOR, 1.0, 1.0, 0.0, ]
  # passed through leading edge of sphere
  in1_obj = [cgo.BEGIN, cgo.LINES, cgo.COLOR, 1.0, 0.0, 0.0, ]
  #in1_obj = [cgo.BEGIN, cgo.POINTS, cgo.COLOR, 1.0, 0.0, 0.0, ]
  # passed through trailing edge of sphere
  in2_obj = [cgo.BEGIN, cgo.LINES, cgo.COLOR, 0.0, 1.0, 0.0, ]
  #in2_obj = [cgo.BEGIN, cgo.POINTS, cgo.COLOR, 0.0, 1.0, 0.0, ]
  asymm_obj = [cgo.BEGIN, cgo.LINES, cgo.COLOR, 1.0, 0.5, 0.5, ]

  for h in range(-cell.hmax,cell.hmax):
    for k in range(-cell.kmax,cell.kmax):
      for l in range(-cell.lmax,cell.lmax):
        (x,y,z) = cell.u_cell_recip.orthogonalize((h,k,l))
        # expand by 10 for ease of viewing in PyMOL
        (x,y,z) = (x*10.,y*10.,z*10.)
        xsq = x*x
        ysq = y*y
        zsq = z*z

        string = [cgo.VERTEX, x,y,z]

        dsq = xsq + ysq + zsq
        
        if dsq <= cell.dminsq_recip:
          if cell.asu.is_inside((h,k,l)):
            asymm_obj.extend(make_cross(cell.r,x,y,z))

          if within_sph1(cell,x,y,z,xsq,ysq,zsq,dsq):
            in1_obj.extend(make_cross(cell.r,x,y,z))
#            in1_obj.extend(string)

          elif within_sph2(cell,x,y,z,xsq,ysq,zsq,dsq):
            in2_obj.extend(make_cross(cell.r,x,y,z))
#            in2_obj.extend(string)

          else:
            in0_obj.extend(string)
        else:
          out_obj.extend(string)


  out_obj.extend([cgo.END,])
  in0_obj.extend([cgo.END,])
  in1_obj.extend([cgo.END,])
  in2_obj.extend([cgo.END,])
  asymm_obj.extend([cgo.END,])
  
  cmd.load_cgo(out_obj,'out')
  cmd.load_cgo(in0_obj,'in0')
  cmd.load_cgo(in1_obj,'in1')
  cmd.load_cgo(in2_obj,'in2')
  cmd.load_cgo(asymm_obj,'asymm')

#if phisweep > 180 - 2 * N.asin(wavelength/(2*dmin)
#draw_axes()
def ewald(cell_param=(30.,40.,80.,90.,99.,90.),sg='p21',dmin=3.0,wavelength=1.5418,phi0=0,phisweep=-90):
  """

  usage:  from within PyMOL, do: 'run ewald.py' to load the program 
  
  and then:

          ewald() to draw the default or specify the options like:
          
          cell_param=(30.,40.,80.,90.,99.,90.)
          sg='p21'
          dmin=3.0
          wavelength=1.5418
                 
  e.g.   ewald(cell_param=(120.,120.,120.,90.,90.,90.),sg='p43',dmin=6.,wavelength=0.99)

  """
  
  Cell = Cell_info(cell_param,sg,dmin,wavelength,phi0,phisweep)
  draw_axes_cyl(Cell)
  draw_ewald0(Cell)
  draw_ewald1(Cell)
  draw_limsph(Cell)
  draw_refl(Cell)

cmd.extend("ewald",ewald)

if __name__ == "__main__":
  import sys

  cell_param = map(float,sys.argv[1:7])
  sg = sys.argv[7]
  dmin = float(sys.argv[8])
  ewald(cell_param,sg,dmin)
