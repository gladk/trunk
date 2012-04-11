# encoding: utf-8
from yade import utils,pack,geom,qt
from pylab import *
import math
qt.View()

O.dt=5e-07
young=2e4


O.materials.append(CohFrictMat(young=8e5,poisson=0.3,density=4e3,frictionAngle=radians(30),normalCohesion=1e5,shearCohesion=1e5,momentRotationLaw=True,label='cylindermat'))

O.materials.append(CohFrictMat(young=8e5,poisson=0.3,density=4e3,frictionAngle=radians(30),normalCohesion=1e5,shearCohesion=1e5,momentRotationLaw=False,label='spheremat'))

rCyl=0.006

cylindersBodies=[]
for i in arange(0,0.3,0.3/15.):
	cylindersBodies.append(utils.chainedCylinder(begin=Vector3(i,0,0),radius=rCyl,end=Vector3(i+0.3/15.,0,0),fixed=False,wire=False,material='cylindermat'))
ChainedState.currentChain=ChainedState.currentChain+1
O.bodies[0].state.blockedDOFs='xyzXYZ'
O.bodies[-1].state.blockedDOFs='xyzXYZ'

O.bodies.append(utils.sphere([0.15,0,2.*rCyl],rCyl,wire=False,fixed=False,material='spheremat'))

def main():
	global Fn,Ft
	IdSphere=len(O.bodies)-1
	if O.iter>50000 :
		O.bodies[IdSphere].dynamic=False
		O.bodies[IdSphere].state.vel[2]=0.1


O.engines=[
	ForceResetter(),
	InsertionSortCollider([
		Bo1_ChainedCylinder_Aabb(),
		Bo1_Sphere_Aabb(),
	]),
	InteractionLoop(
		#Geometric interactions :
		[Ig2_ChainedCylinder_ChainedCylinder_ScGeom6D(),
		Ig2_Sphere_ChainedCylinder_CylScGeom6D(),	#used for the cohesive sphere-cylinder interaction
		Ig2_Sphere_Sphere_ScGeom6D()],
		#Physical interactions
		[Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(setCohesionNow=True,setCohesionOnNewContacts=True,label='ipf'),
		Ip2_FrictMat_FrictMat_FrictPhys()],
		#Lois de contact :
		[Law2_ScGeom6D_CohFrictPhys_CohesionMoment(),
		Law2_CylScGeom6D_CohFrictPhys_CohesionMoment(),	#used for the cohesive sphere-cylinder interaction
		Law2_ScGeom_FrictPhys_CundallStrack(),]
	),
	## Motion equation
	NewtonIntegrator(gravity=(0,0,0),damping=0.3,label='newton'),
	PyRunner(command='main()',iterPeriod=1000),
]



