#!/usr/bin/env python
# encoding: utf-8

# The script implements the show case in this article [Mani2013]
from yade import utils, plot
o = Omega()
fr = 0.5;rho=2000
tc = 0.001; en = 0.7; et = 0.7; 
o.dt = 1.0


r1 = 1.0
r2 = 1.0
Gamma = 20.6*1e-3
Theta = 0
VB3 = 74.2*1e-12


CapillarType = "Lambert"

mat1 = O.materials.append(ViscElCapMat(frictionAngle=fr,density=rho,Vb=VB3,gamma=Gamma,theta=Theta,Capillar=True,CapillarType=CapillarType,tc=tc,en=en,et=et))

d = 0.9999

id1 = O.bodies.append(sphere(center=[0,0,0],                    radius=r2,material=mat1,fixed=True, color=[0,1,0]))
id2 = O.bodies.append(sphere(center=[0,(r1+r2)*d,0],            radius=r2,material=mat1,fixed=True, color=[0,1,0]))
id3 = O.bodies.append(sphere(center=[(r1+r2)*d,(r1+r2)*d,0],    radius=r2,material=mat1,fixed=True, color=[0,1,0]))
id4 = O.bodies.append(sphere(center=[(r1+r2)*d,(r1+r2)*d*2,0],  radius=r2,material=mat1,fixed=True, color=[0,1,0]))
id5 = O.bodies.append(sphere(center=[(r1+r2)*d*2,(r1+r2)*d,0],  radius=r2,material=mat1,fixed=True, color=[0,1,0]))


Vf = 0.5e-1
Vfmin = 0.1e-1

O.bodies[id1].Vf = Vf
O.bodies[id1].Vmin = Vfmin

O.bodies[id2].Vf = Vf
O.bodies[id2].Vmin = Vfmin

O.bodies[id3].Vf = Vf
O.bodies[id3].Vmin = Vfmin

O.bodies[id4].Vf = Vf
O.bodies[id4].Vmin = Vfmin

O.bodies[id5].Vf = Vf
O.bodies[id5].Vmin = Vfmin

vel = 0.0
O.bodies[id1].state.vel=[0,0,vel]
O.bodies[id2].state.vel=[0,0,-vel]
O.bodies[id3].state.vel=[vel,0,0]
O.bodies[id4].state.vel=[-vel,0,0]

o.engines = [
  ForceResetter(),
  InsertionSortCollider([Bo1_Sphere_Aabb()],verletDist=(r1+r2)*5.0),
  InteractionLoop(
    [Ig2_Sphere_Sphere_ScGeom()],
    [Ip2_ViscElCapMat_ViscElCapMat_ViscElCapPhys()],
    [Law2_ScGeom_ViscElCapPhys_Basic()],
  ),
  LiqControl(),
  NewtonIntegrator(damping=0,gravity=[0,0,0]),
  PyRunner(command='showData()',iterPeriod=1),
]


def showData():
  print "Step %d"%O.iter
  print "idB=%d, Vf=%s, Vmin=%s;"%(id1, O.bodies[id1].Vf, O.bodies[id1].Vmin)
  print "idB=%d, Vf=%s, Vmin=%s;"%(id2, O.bodies[id2].Vf, O.bodies[id2].Vmin)
  print "idB=%d, Vf=%s, Vmin=%s;"%(id3, O.bodies[id3].Vf, O.bodies[id3].Vmin)
  print "idB=%d, Vf=%s, Vmin=%s;"%(id4, O.bodies[id4].Vf, O.bodies[id4].Vmin)
  print "idB=%d, Vf=%s, Vmin=%s;"%(id5, O.bodies[id5].Vf, O.bodies[id5].Vmin)
  
  try:
    print "Interaction[1, 2].Vb=%s"%(O.interactions[id1,id2].phys.Vb)
  except:
    pass
  
  try:
    print "Interaction[2, 3].Vb=%s"%(O.interactions[id2,id3].phys.Vb)
  except:
    pass
  
  try:
    print "Interaction[3, 4].Vb=%s"%(O.interactions[id3,id4].phys.Vb)
  except:
    pass
  
  try:
    print "Interaction[3, 5].Vb=%s"%(O.interactions[id3,id5].phys.Vb)
  except:
    pass
  print 


showData()

O.run(1, True)

for i in range(5):
  O.bodies[i].Vf = 0
  O.bodies[i].Vmin = 0

O.interactions[id1,id2].phys.Vb = 1.0
O.interactions[id1,id2].phys.Vmax = 5.0
O.interactions[id2,id3].phys.Vb = 1.0
O.interactions[id2,id3].phys.Vmax = 5.0
O.interactions[id3,id4].phys.Vb = 1.0
O.interactions[id3,id4].phys.Vmax = 5.0
O.interactions[id3,id5].phys.Vb = 1.0
O.interactions[id3,id5].phys.Vmax = 5.0

O.run(1, True)


vel = 1.0
O.bodies[id3].state.vel=[vel,0,0]
O.bodies[id4].state.vel=[vel,0,0]
O.bodies[id5].state.vel=[vel,0,0]

from yade import qt
qt.View()
