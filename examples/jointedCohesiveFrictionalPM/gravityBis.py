# encoding: utf-8

# example of use JCFpm classes : Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM...
# a cubic rock massif, containing a rock joint with ~ 31 deg. dip angle. At one time, the mechanical properties of the joint are degraded, triggering a smooth sliding


# definition of a predicate for use of randomDensePack() function

from yade import pack
dimModele = 10.0
pred = pack.inAlignedBox((0,0,0),(dimModele,dimModele,dimModele))


# the packing of spheres :

def mat(): return JCFpmMat(type=1,young=1e8,poisson=0.3,frictionAngle=radians(30),density=3000,tensileStrength=1e6,cohesion=1e6,jointNormalStiffness=1e7,jointShearStiffness=1e7,jointCohesion=1e6,jointFrictionAngle=radians(20),jointDilationAngle=0.0)
nSpheres = 3000.0
poros=0.13
rMeanSpheres = dimModele * pow(3.0/4.0*(1-poros)/(pi*nSpheres),1.0/3.0)
print ''
print 'Creating a cubic sample of spheres (may take some time and cause warnings)'
print ''
sp = pack.randomDensePack(pred,radius=rMeanSpheres,rRelFuzz=0.3,memoizeDb='/tmp/gts-triax-packings.sqlite',returnSpherePack=True)
sp.toSimulation(color=(0.9,0.8,0.6),wire=False,material=mat)
print ''
print 'Sample created !'


# Definition of the facets for joint's geometry

import gts
# joint with ~ 31 deg. dip angle
v1 = gts.Vertex(0 , 0 , 0.8*dimModele)
v2 = gts.Vertex(0 , dimModele , 0.8*dimModele )
v3 = gts.Vertex(dimModele , dimModele , 0.2*dimModele)
v4 = gts.Vertex(dimModele , 0 , 0.2*dimModele)

e1 = gts.Edge(v1,v2)
e2 = gts.Edge(v2,v4)
e3 = gts.Edge(v4,v1)
f1 = gts.Face(e1,e2,e3)

e4 = gts.Edge(v4,v3)
e5 = gts.Edge(v3,v2)
f2 = gts.Face(e2,e4,e5)

s1 = gts.Surface()
s1.add(f1)
s1.add(f2)

facet = gtsSurface2Facets(s1,wire = False,material=mat)
O.bodies.append(facet)

yade.qt.View()
yade.qt.Controller()

O.saveTmp()
# identification of spheres onJoint, and so on:
execfile('identifBis.py')
dim=utils.aabbExtrema()
xsup=dim[1][0]
yinf=dim[0][1]
ysup=dim[1][1]
zinf=dim[0][2]
zsup=dim[1][2]

e=2*rMeanSpheres

for o in O.bodies:
   if isinstance(o.shape,Sphere):
      o.shape.color=(0.9,0.8,0.6)
      ## to fix boundary particles on ground
      if o.state.pos[2]<(zinf+2*e) :
	 o.state.blockedDOFs+='xyz'
	 o.shape.color=(1,1,1)

      ## to identify indicator on top
      if o.state.pos[2]>(zsup-e) and o.state.pos[0]>(xsup-e) and o.state.pos[1]>((yinf+ysup-e)/2.0) and o.state.pos[1]<((yinf+ysup+e)/2) : 
	refPoint=o.id

O.bodies[refPoint].shape.highlight=True

#### Engines definition
O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(),]),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom()],
#		[Ip2_JCFpmMat_JCFpmMat_JCFpmPhys(cohesiveTresholdIteration=1,alpha=0.3,tensileStrength=1e6,cohesion=1e6,jointNormalStiffness=1e7,jointShearStiffness=1e7,jointCohesion=1e6,jointFrictionAngle=radians(20),jointDilationAngle=0.0)],
		[Ip2_JCFpmMat_JCFpmMat_JCFpmPhys(cohesiveTresholdIteration=1)],
		[Law2_ScGeom_JCFpmPhys_JointedCohesiveFrictionalPM(smoothJoint=True)]
	),
	GlobalStiffnessTimeStepper(timestepSafetyCoefficient=0.8),
	PyRunner(iterPeriod=1000,initRun=False,command='jointStrengthDegradation()'),
	PyRunner(iterPeriod=10,initRun=True,command='dataCollector()'),
	NewtonIntegrator(damping=0.7,gravity=(0.,0.,-10.)),

]

#### dataCollector
from yade import plot
plot.plots={'iterations':'v','x':'z'}

def dataCollector():
	R=O.bodies[refPoint]
	plot.addData(v=R.state.vel[2],z=R.state.pos[2],x=R.state.pos[0],iterations=O.iter,t=O.realtime)

#### joint strength degradation
stableIter=1000
stableVel=0.001
degrade=True
def jointStrengthDegradation():
    global degrade
    if degrade and O.iter>=stableIter and abs(O.bodies[refPoint].state.vel[2])<stableVel :
	print 'Equilibrium reached \nJoint cohesion canceled now !', ' | iteration=', O.iter
	degrade=False
	for i in O.interactions:
	    if i.phys.isOnJoint : 
		if i.phys.isCohesive:
		  i.phys.isCohesive=False
		  i.phys.FnMax=0.
		  i.phys.FsMax=0.

print 'Seeking after an initial equilibrium state'
print ''
O.run(10000)
plot.plot()# note the straight trajectory (z(x) plot)during sliding step (before free fall) despite the discretization of joint plane with spheres
