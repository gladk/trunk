"""
This example shows how LSMGenGeo generates Cylinder 
with given parameters.

Be sure LSMGenGeo library is installed.

The result is:
	2 files:
		"cyl.geo" is the geometry file which can be imported into YADE with ymport.gengeoFile() function
		"cyl.vtk" is the VTK-filed which can be opened by any VTK-based software, for example Paraview
	spheres, imported into the YADE simulation, according to generated geometry

http://www.access.edu.au/lsmgengeo_python_doc/current/pythonapi/html/GenGeo-module.html
https://svn.esscc.uq.edu.au/svn/esys3/lsm/contrib/LSMGenGeo/
"""
from yade import pack,ymport

try:
	from GenGeo import *
except ImportError:
	from gengeo import *

import sys

fileName="cyl"
radius=15
length=100

minRadius=0.5
maxRadius=2.0

origin=Vector3(0.0,0.0,0.0)
axis=Vector3(0.0,0.0,1.0)

minPoint=Vector3(-radius,-radius,-length)
maxPoint=Vector3(radius,radius,length)

box=CylinderVol (
  origin = origin,
  axis = axis,
  radius = radius,
  length = 2.0*length
)

mntable=MNTable3D (
  minPoint=minPoint,
  maxPoint=maxPoint,
  gridSize=2.5*maxRadius,
  numGroups=1
)

packer=InsertGenerator3D(
  minRadius=minRadius,
  maxRadius=maxRadius,
  insertFails=1000,
  maxIterations=1000,
  tolerance=1.0e-6,
  seed=0
)

packer.generatePacking (
  volume=box,
  ntable=mntable,
  groupID=0
)

#Write data to files
mntable.write(fileName+".geo",1)
mntable.write(fileName+".vtu",2)

#Add material
O.materials.append(FrictMat(young=1e9,poisson=.25,frictionAngle=0.5,density=1e3))

#Parameters, which will be passed into spheres and facets creators
kw={'material':0}

#Import the GenGeo geometry directly into the YADE simulation
O.bodies.append(ymport.gengeo(mntable,shift=[-1.0,-1.0,-1.0],scale=2.0,color=(1,1,0),**kw))

try:
	from yade import qt
	qt.Controller()
	qt.View()
except ImportError: pass

