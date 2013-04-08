# encoding: utf-8
# 2010 Chiara Modenese <c.modenese@gmail.com>

# Script to test the contact damping in HM (both in the normal and shear direction)



#__________________________________________________________________
# Geometry sphere and box
r2=1e-2 # radii [m]
p2=[r2,0,0] # center positions [m]
center=[-r2/2.,0,0] # center [m]
extents=[r2/2.,3*r2,3*r2] # half edge lenght [m]

#__________________________________________________________________
# Material
young=600.0e6 # [N/m^2]
poisson=0.6 
density=2.60e3 # [kg/m^3]
frictionAngle=26 # [°]

# Append geometry and material
O.materials.append(FrictMat(young=young,poisson=poisson,density=density,frictionAngle=frictionAngle))
O.bodies.append(box(center=center,extents=extents,fixed=True,wire=True)) # body id=0
O.bodies.append(sphere(p2,r2,fixed=True,wire=True)) # body id=1
m_sphere=O.bodies[1].state.mass
O.bodies[0].state.mass=m_sphere # set the mass of the box the same as the mass of the sphere
O.bodies[1].state.blockedDOFs='XYZ' # block particles rotations

#__________________________________________________________________
# list of engines
O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb()]),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
		[Ip2_FrictMat_FrictMat_MindlinPhys(label='damping')],
		[Law2_ScGeom_MindlinPhys_Mindlin(label='contactLaw')]
	),
	ForceEngine(force=(-20,0,0),ids=[1],label='force'),	
	NewtonIntegrator(damping=0.0),
	PyRunner(iterPeriod=1,command='myAddPlotData()'),
]

#__________________________________________________________________
# define contact damping coefficients
damping.betan=0.5 # normal direction
damping.betas=0.0 # shear direction
damping.useDamping=True

#__________________________________________________________________
# time step
O.dt=.2*PWaveTimeStep()
O.saveTmp('init')

#__________________________________________________________________
from yade import qt
qt.View()
qt.Controller()

#__________________________________________________________________
# plot some results
from yade import plot

plot.labels=dict(Fn='$Normal\,force$',un='$Overlapping$',time='$Time$',time_='$Time$',Fs='$Shear\,force$',t='$Time$')
plot.plots={'time':('un',),'time_':('Fn',),'t':('Fs')}

def myAddPlotData():
	i=O.interactions[0,1]
	plot.addData(Fn=i.phys.normalForce[0],Fs=i.phys.shearForce[1],un=i.geom.penetrationDepth,time=O.time,time_=O.time,t=O.time)

# We will have:
# 1) data in graphs (if you call plot.plot())
# 2) data in file (if you call plot.saveGnuplot('/home/chia/Desktop/Output/out')
# 3) data in memory as plot.data['un'], plot.data['fn'], etc. under the labels they were saved

O.run(1000,True)

# Now test the shear direction, having previously obtained a constant normal overlapping 
damping.betas=0.5
force.force=(-20,5,0) # assign a force also in the shear direction
O.run(1000,True)
plot.plot(subPlots=False)

# NOTE about the results:
# In the graphs, you will see the behaviour both of the normal and shear force. As expected from the damped solution, their values oscillate around the equilibrium position, which is eventually reached after few iterations. 




