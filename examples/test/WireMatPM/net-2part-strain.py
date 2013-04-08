# -*- coding: utf-8 -*-
# encoding: utf-8
from yade import  ymport, qt


## definition of some colors for colored text output in terminal
BLUE = '\033[94m'
GREEN = '\033[92m'
YELLOW = '\033[93m'
RED = '\033[91m'
BLACK = '\033[0m'

#### short description of script
print BLUE+'Simple test for two particles to test contact law with '+RED+'UniaxialStrainer'+BLUE+'.'+BLACK

#### define parameters for the net
# mesh opening size
mos = 80./1000.
a = mos/sqrt(3)
# wire diameter
d = 2.7/1000.
# particle radius
radius = d*5.
# define piecewise lineare stress-strain curve
strainStressValues=[(0.0019230769,2.5e8),(0.0192,3.2195e8),(0.05,3.8292e8),(0.15,5.1219e8),(0.25,5.5854e8),(0.3,5.6585e8),(0.35,5.6585e8)]
# elastic material properties
particleVolume = 4./3.*pow(radius,3)*pi
particleMass = 3.9/1000.
density = particleMass/particleVolume
young = strainStressValues[0][1] / strainStressValues[0][0]
poisson = 0.3


#### material definition
netMat = O.materials.append(WireMat(young=young,poisson=poisson,frictionAngle=radians(30),density=density,isDoubleTwist=False,diameter=d,strainStressValues=strainStressValues,lambdaEps=0.4,lambdak=0.21))


#### create boddies, default: dynamic=True
O.bodies.append( sphere([0,0,0], radius, wire=False, color=[1,0,0], highlight=False, material=netMat) )
O.bodies.append( sphere([0,a,0], radius, wire=False, color=[0,1,0], highlight=False, material=netMat) )

FixedSphere=O.bodies[0]
MovingSphere=O.bodies[1]

FixedSphere.dynamic=True
MovingSphere.dynamic=True


#### initialize values for UniaxialStrainer
bb = uniaxialTestFeatures(axis=1)
negIds,posIds,axis,crossSectionArea=bb['negIds'],bb['posIds'],bb['axis'],bb['area']
strainRateTension = 1./a
setSpeeds = True


#### define simulation to create link
interactionRadius=2.
O.engines = [
	ForceResetter(),
	InsertionSortCollider( [Bo1_Sphere_Aabb(aabbEnlargeFactor=interactionRadius,label='aabb')] ), 

	InteractionLoop(
	[Ig2_Sphere_Sphere_ScGeom(interactionDetectionFactor=interactionRadius,label='Ig2ssGeom')],
	[Ip2_WireMat_WireMat_WirePhys(linkThresholdIteration=1,label='interactionPhys')],
	[Law2_ScGeom_WirePhys_WirePM(linkThresholdIteration=1,label='interactionLaw')]
	),
	NewtonIntegrator(damping=0.)
]


#### create link (no time step needed since loading is involved in this step)
O.step() # create cohesive link (cohesiveTresholdIteration=1)


#### initializes now the interaction detection factor
aabb.aabbEnlargeFactor=-1.
Ig2ssGeom.interactionDetectionFactor=-1.

## time step definition
O.dt = 1e-5
## critical time step proposed by Bertrand
#O.dt = 0.2*sqrt(particleMass/(2.*O.interactions[0,1].phys.kn))

#### plot some results
from yade import plot

plot.plots={'un':('Fn',)}
plot.plot(noShow=False, subPlots=False)

def addPlotData():
	try:
		i=O.interactions[FixedSphere.id,MovingSphere.id]
		plot.addData( Fn=i.phys.normalForce.norm(), un=(O.bodies[1].state.pos[1]-O.bodies[0].state.pos[1])-a )
		#plot.saveGnuplot('net-2part-strain')
	except:
		print "No interaction!"
		O.pause()


#### define simulation
O.engines += [UniaxialStrainer(strainRate=strainRateTension,axis=axis,asymmetry=1,posIds=posIds,negIds=negIds,crossSectionArea=crossSectionArea,blockDisplacements=True,blockRotations=True,setSpeeds=setSpeeds,label='strainer')] + [PyRunner(initRun=True,iterPeriod=1,command='addPlotData()')]


#### to see it
v=qt.Controller()
v=qt.View()
rr=qt.Renderer()
rr.intrAllWire=True
