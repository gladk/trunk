#--- bruno.chareyre@hmg.inpg.fr ---
#!/usr/bin/python
# -*- coding: utf-8 -*-

# Experiment beam-like behaviour with chained cylinders + CohFrict connexions

from yade import pack

young=4.0e6
poisson=3
density=1e3
frictionAngle1=radians(15)
frictionAngle2=radians(15)
frictionAngle3=radians(15)
O.dt=5e-4

O.materials.append(FrictMat(young=4000000.0,poisson=0.5,frictionAngle=frictionAngle1,density=1600,label='spheremat'))
O.materials.append(FrictMat(young=1.0e6,poisson=0.2,density=2.60e3,frictionAngle=frictionAngle2,label='walllmat'))

Ns=90
sp=pack.SpherePack()

#cohesive spheres crash because sphere-cylnder functor geneteragets ScGeom3D
#O.materials.append(CohFrictMat(young=1.0e5,poisson=0.03,density=2.60e2,frictionAngle=frictionAngle,label='spheremat'))

if os.path.exists("cloud4cylinders"+`Ns`):
 	print "loading spheres from file"
	sp.load("cloud4cylinders"+`Ns`)
else:
 	print "generating spheres"
 	Ns=sp.makeCloud(Vector3(-0.3,0.2,-1.0),Vector3(+0.3,+0.5,+1.0),-1,.2,Ns,False,0.8)
	sp.save("cloud4cylinders"+`Ns`)

O.bodies.append([sphere(center,rad,material='spheremat') for center,rad in sp])
walls=aabbWalls((Vector3(-0.3,-0.15,-1),Vector3(+0.3,+1.0,+1)),thickness=.1,material='walllmat')
wallIds=O.bodies.append(walls)

O.initializers=[
	BoundDispatcher([Bo1_Sphere_Aabb(),Bo1_ChainedCylinder_Aabb(),Bo1_Box_Aabb()])
]

O.engines=[
	ForceResetter(),
	InsertionSortCollider([
		Bo1_ChainedCylinder_Aabb(),
		Bo1_Sphere_Aabb(),
		Bo1_Box_Aabb()
	]),
	InteractionLoop(
		[Ig2_ChainedCylinder_ChainedCylinder_ScGeom6D(), Ig2_Sphere_ChainedCylinder_CylScGeom(), Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],			
		[Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(setCohesionNow=True,setCohesionOnNewContacts=False,label='ipf'),Ip2_FrictMat_FrictMat_FrictPhys()],
		[Law2_ScGeom6D_CohFrictPhys_CohesionMoment(label='law'),Law2_ScGeom_FrictPhys_CundallStrack(),Law2_CylScGeom_FrictPhys_CundallStrack()]
	),
	## Motion equation
	NewtonIntegrator(damping=0.05,label='newton',gravity=[1,-9.81,0]),
]

#Assemble cylinders in sinusoidal shapes
O.materials.append(CohFrictMat(young=young,poisson=poisson,density=3.0*density,frictionAngle=frictionAngle3,normalCohesion=1e40,shearCohesion=1e40,momentRotationLaw=True,label='cylindermat'))
Ne=30
dy=0.03
dx=0.2
dz=0.2
Nc=1 #nb. of additional chains
for j in range(-Nc, Nc+1):
	dyj = abs(float(j))*dy
	dxj = abs(float(j))*dx
	dzj = float(j)*dz
	for i in range(0, Ne):
		omega=20/float(Ne); hy=0.0; hz=0.05; hx=1.5;
		px=float(i)*hx/float(Ne)-0.8+dxj; py=sin(float(i)*omega)*hy+dyj; pz=cos(float(i)*omega)*hz+dzj;
		px2=float(i+1.)*hx/float(Ne)-0.8+dxj; py2=sin(float(i+1.)*omega)*hy+dyj; pz2=cos(float(i+1.)*omega)*hz+dzj;
		chainedCylinder(begin=Vector3(pz,py,px), radius=0.02,end=Vector3(pz2,py2,px2),color=Vector3(0.6,0.5,0.5),material='cylindermat')
		if (i == Ne-1): #close the chain with a node of size 0
			print "closing chain"
			b=chainedCylinder(begin=Vector3(pz2,py2,px2), radius=0.02,end=Vector3(pz2,py2,px2),color=Vector3(0.6,0.5,0.5),material='cylindermat')
			b.state.blockedDOFs='xyzXYZ'
	ChainedState.currentChain=ChainedState.currentChain+1

O.saveTmp()
