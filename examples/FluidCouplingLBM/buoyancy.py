##############################################################################################################################
#Authors: Luc Sibille  luc.sibille@3sr-grenoble.fr
#         Franck Lomine  franck.lomine@insa-rennes.fr
#
#   Script to simulate buoyancy in a granular assembly with the DEM-LBM coupling.
#   This script was written for a practical work in the OZ ALERT school in Grenoble 2011.
#   Script updated in 2014
##############################################################################################################################


from yade import pack,timing,utils


rMean = 0.00035  #mean radius (m) of solid particles

Gravity=-0.1 #gravity (m.s-2) that will be applied to solid particles

#definition of the simulation domain (position of walls), the domain is a 3D domain but particles will be generated on a unique plane at z=0, be carefull that wall normal to z directions do not touch particles.
lowerCornerW = (0.00001,0.00001,-0.002)
upperCornerW = (0.00701,0.00501,0.002)

#definitions of the domain for particles: positions in z direction is set to 0 to force a 2D granular assembly 
lowerCornerS = (0.00001,0.00001,0)
upperCornerS = (0.00701,0.00501,0)


nbSpheres = 70 #this not the real number of particles but number of times where we try to insert a new particle, the real number of particles can be much lower! 



## Build of the engine vector
O.engines=[

	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb()]),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
		[Ip2_FrictMat_FrictMat_FrictPhys()],
		[Law2_ScGeom_FrictPhys_CundallStrack()]
	),

	
	HydrodynamicsLawLBM(	EngineIsActivated=False,
				WallYm_id=0,
				WallYp_id=1,
				WallXm_id=2,
				WallXp_id=3, 
				WallZp_id=5, 
				WallZm_id=4,
				useWallYm=1,
				useWallYp=1,
				YmBCType=2,
				YpBCType=2,

				useWallXm=0,
				useWallXp=0,
				XmBCType=1,
				XpBCType=1,
				LBMSavedData='spheres,velXY,rho,nodeBD',
				tau=1.1,
				dP=(0,0,0),
				IterSave=200,
				IterPrint=100,
				RadFactor=0.6,  # The radius of particles seen by the LBM engine is reduced here by a factor RadFactor=0.6 to allow flow between particles in this 2D case.
				Nx=250,
				Rho=1000,
				Nu=1.0e-6,
				periodicity='',
				bc='',
				VbCutOff=0.0,
				applyForcesAndTorques=True,
				label="Elbm"
	),


	NewtonIntegrator(damping=0.2,gravity=(Gravity,0.,0.),label="ENewton"),
]



###############################################################################################################
# Creation of 6 walls at the limit of the simulation domain

# defintiion of the material for walls
O.materials.append(FrictMat(young=50e6,poisson=.5,frictionAngle=0.0,density=3000,label='walls'))

## create walls around the packing
wallOversizeFactor=1.001
thickness=0.00001;
#bottom box
center= ((lowerCornerW[0]+upperCornerW[0])/2,lowerCornerW[1]-thickness/2.0,(lowerCornerW[2]+upperCornerW[2])/2)
halfSize= (wallOversizeFactor*fabs(lowerCornerW[0]-upperCornerW[0])/2+thickness,thickness/2.0,wallOversizeFactor*fabs(lowerCornerW[2]-upperCornerW[2])/2+thickness)
b1=utils.box(center=[center[0],center[1],center[2]],extents=[halfSize[0],halfSize[1],halfSize[2]],color=[0,1,0],fixed=True,wire=True,material='walls')
O.bodies.append(b1)
#--
#Top box
center=((lowerCornerW[0]+upperCornerW[0])/2,upperCornerW[1]+thickness/2.0,(lowerCornerW[2]+upperCornerW[2])/2)
halfSize =(wallOversizeFactor*fabs(lowerCornerW[0]-upperCornerW[0])/2+thickness,thickness/2.0,wallOversizeFactor*fabs(lowerCornerW[2]-upperCornerW[2])/2+thickness)
b2=utils.box(center=[center[0],center[1],center[2]],extents=[halfSize[0],halfSize[1],halfSize[2]],color=[0,1,0],fixed=True,wire=True,material='walls')
O.bodies.append(b2)
#--
center=(lowerCornerW[0]-thickness/2.0,(lowerCornerW[1]+upperCornerW[1])/2,(lowerCornerW[2]+upperCornerW[2])/2)
halfSize=(thickness/2.0,wallOversizeFactor*fabs(lowerCornerW[1]-upperCornerW[1])/2+thickness,wallOversizeFactor*fabs(lowerCornerW[2]-upperCornerW[2])/2+thickness)
b3=utils.box(center=[center[0],center[1],center[2]],extents=[halfSize[0],halfSize[1],halfSize[2]],color=[0,1,0],fixed=True,wire=True,material='walls')
O.bodies.append(b3)
#--
center=(upperCornerW[0]+thickness/2.0,(lowerCornerW[1]+upperCornerW[1])/2,(lowerCornerW[2]+upperCornerW[2])/2)
halfSize=(thickness/2.0,wallOversizeFactor*fabs(lowerCornerW[1]-upperCornerW[1])/2+thickness,wallOversizeFactor*fabs(lowerCornerW[2]-upperCornerW[2])/2+thickness)
b4=utils.box(center=[center[0],center[1],center[2]],extents=[halfSize[0],halfSize[1],halfSize[2]],color=[0,1,0],fixed=True,wire=True,material='walls')
O.bodies.append(b4)
#--
center=((lowerCornerW[0]+upperCornerW[0])/2,(lowerCornerW[1]+upperCornerW[1])/2,lowerCornerW[2]-thickness/2.0)
halfSize=(wallOversizeFactor*fabs(lowerCornerW[0]-upperCornerW[0])/2+thickness,wallOversizeFactor*fabs(lowerCornerW[1]-upperCornerW[1])/2+thickness,thickness/2.0)
b5=utils.box(center=[center[0],center[1],center[2]],extents=[halfSize[0],halfSize[1],halfSize[2]],color=[0,1,0],fixed=True,wire=True,material='walls')
O.bodies.append(b5)
#--
center=((lowerCornerW[0]+upperCornerW[0])/2,(lowerCornerW[1]+upperCornerW[1])/2,upperCornerW[2]+thickness/2.0)
halfSize=(wallOversizeFactor*fabs(lowerCornerW[0]-upperCornerW[0])/2+thickness,wallOversizeFactor*fabs(lowerCornerW[1]-upperCornerW[1])/2+thickness,thickness/2.0);	 	
b6=utils.box(center=[center[0],center[1],center[2]],extents=[halfSize[0],halfSize[1],halfSize[2]],color=[0,1,0],fixed=True,wire=True,material='walls')
O.bodies.append(b6)


###############################################################################################################
# Creation of the assembly of particles

# defintiion of the material for particles
O.materials.append(FrictMat(young=50e6,poisson=.5,frictionAngle=35*3.1415926535/180,density=3000,label='spheres'))
## use a SpherePack object to generate a random loose particles packing
sp=pack.SpherePack()
sp.makeCloud(lowerCornerS,upperCornerS,-1,0.33,nbSpheres,False, 0.5,seed=1) #"seed" make the "random" generation always the same
sp.toSimulation(material='spheres')


# loop over bodies to change their 3D inertia into 2D inertia (by default in YADE particles are 3D spheres, here we want to cylinder with a length=1).
for s in O.bodies:
	if isinstance(s.shape,Box): continue	

	r=s.shape.radius
	oldm=s.state.mass
	oldI=s.state.inertia

	m=oldm*3./4./r
	s.state.mass=m

	s.state.inertia[0] = 15./16./r*oldI[0]	#inertia with respect to x and y axes are not used and the computation here is wrong
	s.state.inertia[1] = 15./16./r*oldI[1]  #inertia with respect to x and y axes are not used and the computation here is wrong
	s.state.inertia[2] = 15./16./r*oldI[2]  #only inertia with respect to z axis is usefull



O.dt=0.00005  #use a fix value of time step, it is better to not use a time stepper computing a new step at each iteration since DEM time step is eventually fixed by the LBM engine.

yade.qt.View() #to see what is happening 

print"___________________"
print"PHASE 1"
print "Settlement of particle under gravity only, fluid flow is not activated at the present time... please wait"
print"___________________"

O.run(60000,1)  #settlement of particle under gravity only, no action of the fluid flow.



#definition of a runnable python function to increase progressively the fluid pressure gradient.
def IncrDP():
	currentDP=Elbm.dP
	Elbm.dP=(currentDP[0]+1.0/10000.0,0,0)
	#O.engines[6].dP=currentDP+1.0/100.0	
	if not(O.iter%100):
		print "dP=", Elbm.dP[0]

O.engines=O.engines+[PyRunner(iterPeriod=1,command='IncrDP()')]

#Necessary to initiate correctly the LBM engine
O.resetTime()

#Activation of the LBM engine
Elbm.EngineIsActivated=True

#Cundall damping for solid particles removed (the fluid viscosity should be enough to damp the solid particles ... in general, not always...)
ENewton.damping=0.0

# Now you just need to run the simulation as long as you want to ... note that if you wait to long the pressure gradient will become very large inducing great fluid velocities, and the fluid simulation may become crazy!! Look at the Mach number (M) printed in the console, it should not be too high...

print"___________________"
print"PHASE 2"
print "Fluid flow has been activated now with a gradual increase of the pressure gradient."
print "Run the simulation, until the pressure gradient is strong enough to compensate the gravity forces, then particles will mouve to the right."
print"..."
print"..."
print "Velocity and pressure field of the fluid can be plotted with the matlab script yadeLBMvisu.m"
print "Terzaghi critical hydraulic gradient can be compared with the one deduced from the simulation with the matlab script Buoyancy_Gradient.m"
print"___________________"


