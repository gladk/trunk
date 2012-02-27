#!/usr/bin/python
# -*- coding: utf-8 -*-

from yade import utils
from yade import ymport

## Omega
o=Omega() 

## PhysicalParameters 
Density=2400
frictionAngle=radians(35)
sphereRadius=0.05
tc = 0.001
en = 0.3
es = 0.3

## Import wall's geometry
params=utils.getViscoelasticFromSpheresInteraction(tc,en,es)
facetMat=O.materials.append(ViscElMat(frictionAngle=frictionAngle,**params)) # **params sets kn, cn, ks, cs
sphereMat=O.materials.append(ViscElMat(density=Density,frictionAngle=frictionAngle,**params))

walls = O.bodies.append(ymport.stl('ring.stl',material=facetMat))

def fill_cylinder_with_spheres(sphereRadius,cylinderRadius,cylinderHeight,cylinderOrigin,cylinderSlope):
	spheresCount=0
	for h in xrange(0,cylinderHeight/sphereRadius/2):
			for r in xrange(1,cylinderRadius/sphereRadius/2):
				dfi = asin(0.5/r)*2
				for a in xrange(0,6.28/dfi):
					x = cylinderOrigin[0]+2*r*sphereRadius*cos(dfi*a)
					y = cylinderOrigin[1]+2*r*sphereRadius*sin(dfi*a)
					z = cylinderOrigin[2]+h*2*sphereRadius
					s=utils.sphere([x,y*cos(cylinderSlope)+z*sin(cylinderSlope),z*cos(cylinderSlope)-y*sin(cylinderSlope)],sphereRadius,material=sphereMat)
					o.bodies.append(s)
					spheresCount+=1
	return spheresCount

# Spheres
spheresCount=0
spheresCount+=fill_cylinder_with_spheres(sphereRadius,0.5,0.10,[0,0,0],radians(0))
print "Number of spheres: %d" % spheresCount


## Engines 
o.engines=[
	## Resets forces and momenta the act on bodies
	ForceResetter(),

	## Using bounding boxes find possible body collisions.
	InsertionSortCollider([
		Bo1_Sphere_Aabb(),
		Bo1_Facet_Aabb(),
	]),
	# Interactions
	InteractionLoop(
		## Create geometry information about each potential collision.
		[Ig2_Sphere_Sphere_ScGeom(), Ig2_Facet_Sphere_ScGeom()],
		## Create physical information about the interaction.
		[Ip2_ViscElMat_ViscElMat_ViscElPhys()],
		## Constitutive law
		[Law2_ScGeom_ViscElPhys_Basic()],
	),
	## Apply gravity
	GravityEngine(gravity=[0,-9.81,0]),
	## Cundall damping must been disabled!
	NewtonIntegrator(damping=0),
	## Apply kinematics to walls
    ## angularVelocity = 0.73 rad/sec = 7 rpm
	RotationEngine(subscribedBodies=walls,rotationAxis=[0,0,1],rotateAroundZero=True,angularVelocity=0.73)

]

# FIXME: why this not worked properly?
#for b in O.bodies:
    #if isinstance(b.shape.name,Sphere):
		 #b.state.blockedDOFs=['z']

for b in o.bodies:
	if b.shape.name=='Sphere':
		b.state.blockedDOFs=['z']


O.dt=0.02*tc

O.saveTmp('init');

from yade import qt
renderer=qt.Renderer()
renderer.wire=True
#qt.Controller()
qt.View()
O.run()

