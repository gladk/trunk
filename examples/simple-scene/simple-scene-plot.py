#!/usr/bin/python
# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('TkAgg')

O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb()]),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
		[Ip2_FrictMat_FrictMat_FrictPhys()],
		[Law2_ScGeom_FrictPhys_CundallStrack()]
	),
	NewtonIntegrator(damping=.2,gravity=(0,0,-9.81)),
	###
	### NOTE this extra engine:
	###
	### You want snapshot to be taken every 1 sec (realTimeLim) or every 50 iterations (iterLim),
	### whichever comes soones. virtTimeLim attribute is unset, hence virtual time period is not taken into account.
	PyRunner(iterPeriod=20,command='myAddPlotData()')
]
O.bodies.append(box(center=[0,0,0],extents=[.5,.5,.5],fixed=True,color=[1,0,0]))
O.bodies.append(sphere([0,0,2],1,color=[0,1,0]))
O.dt=.002*PWaveTimeStep()


############################################
##### now the part pertaining to plots #####
############################################

from yade import plot
## we will have 2 plots:
## 1. t as function of i (joke test function)
## 2. i as function of t on left y-axis ('|||' makes the separation) and z_sph, v_sph (as green circles connected with line) and z_sph_half again as function of t
plot.plots={'i':('t'),'t':('z_sph',None,('v_sph','go-'),'z_sph_half')}

## this function is called by plotDataCollector
## it should add data with the labels that we will plot
## if a datum is not specified (but exists), it will be NaN and will not be plotted
def myAddPlotData():
	sph=O.bodies[1]
	## store some numbers under some labels
	plot.addData(t=O.time,i=O.iter,z_sph=sph.state.pos[2],z_sph_half=.5*sph.state.pos[2],v_sph=sph.state.vel.norm())
print "Now calling plot.plot() to show the figures. The timestep is artificially low so that you can watch graphs being updated live."
plot.liveInterval=.2
plot.plot(subPlots=False)
O.run(int(2./O.dt));
#plot.saveGnuplot('/tmp/a')
## you can also access the data in plot.data['i'], plot.data['t'] etc, under the labels they were saved.
