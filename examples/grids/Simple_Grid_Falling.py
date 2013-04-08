# encoding: utf-8

from yade import pack,geom,qt
from pylab import *
qt.View()


O.engines=[
	ForceResetter(),
	InsertionSortCollider([
		Bo1_GridConnection_Aabb(),
	]),
	InteractionLoop(
		[Ig2_GridNode_GridNode_GridNodeGeom6D()],
		[Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(setCohesionNow=True,setCohesionOnNewContacts=False)],
		[Law2_ScGeom6D_CohFrictPhys_CohesionMoment()]
	),
	NewtonIntegrator(gravity=(0,0,-10),damping=0.1,label='newton')
]


O.materials.append(CohFrictMat(young=3e2,poisson=0.3,density=1e1,frictionAngle=10,normalCohesion=1e7,shearCohesion=1e7,momentRotationLaw=True,label='spheremat'))


### Parameters of a rectangular grid ###
L=0.1 #length [m]
l=0.05	#width	[m]
nbL=10	#number of nodes for the length	[#]
nbl=5	#number of nodes for the width	[#]
r=L/100.	#radius
color=[255./255.,102./255.,0./255.]
nodesIds=[]
#Create all nodes first :
for i in range(0,nbL):
	for j in range(0,nbl):
		nodesIds.append( O.bodies.append(gridNode([i*L/nbL,j*l/nbl,0],r,wire=False,fixed=False,material='spheremat',color=color)) )

#Create connection between the nodes
for i in range(0,len(nodesIds)):
	for j in range(i+1,len(nodesIds)):
		dist=(O.bodies[i].state.pos - O.bodies[j].state.pos).norm()
		if(dist<=L/nbL*1.01):
			O.bodies.append( gridConnection(i,j,r,color=color) )

#Set a fixed node
O.bodies[0].dynamic=False

O.dt=1e-05
O.saveTmp()


