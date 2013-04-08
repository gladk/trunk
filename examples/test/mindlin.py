#!/usr/bin/python
# -*- coding: utf-8 -*-
##
## SCRIPT TO TEST A NEW CONSTITUTIVE LAW (MINDLIN - nonlinear elastic model)

## list of engines
O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb()]),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(),Ig2_Box_Sphere_ScGeom()],
		[Ip2_FrictMat_FrictMat_MindlinPhys()],
		[Law2_ScGeom_MindlinPhys_Mindlin()]
	),
	NewtonIntegrator(damping=0.0,gravity=(10,0,0)),
	###
	### NOTE this extra engine:
	###
	### You want snapshot to be taken every 1 sec (realTimeLim) or every 50 iterations (iterLim),
	### whichever comes soones. virtTimeLim attribute is unset, hence virtual time period is not taken into account.
	PyRunner(iterPeriod=1,command='myAddPlotData()')
]

## define and append material
mat=FrictMat(young=600.0e6,poisson=0.6,density=2.60e3,frictionAngle=26,label='Friction')
O.materials.append(mat)

## create two spheres (one will be fixed) and append them

s0=sphere([0,0,0],1,color=[0,1,0],fixed=False,wire=True,material='Friction')
s1=sphere([2,0,0],1,color=[0,2,0],fixed=True,wire=True,material='Friction')
O.bodies.append(s0)
O.bodies.append(s1)

## time step
O.dt=.2*PWaveTimeStep()
O.saveTmp('Mindlin')

from yade import qt
qt.View()
qt.Controller()

############################################
##### now the part pertaining to plots #####
############################################

from yade import plot
## make one plot: step as function of fn
plot.plots={'un':('fn')}

## this function is called by plotDataCollector
## it should add data with the labels that we will plot
## if a datum is not specified (but exists), it will be NaN and will not be plotted

def myAddPlotData():
	if O.interactions[0,1].isReal:
		i=O.interactions[0,1]
		## store some numbers under some labels
		plot.addData(fn=i.phys.normalForce[0],step=O.iter,un=2*s0.shape.radius-s1.state.pos[0]+s0.state.pos[0],kn=i.phys.kn)	

O.run(100,True); plot.plot(subPlots=False)

## We will have:
## 1) data in graphs (if you call plot.plot())
## 2) data in file (if you call plot.saveGnuplot('/tmp/a')
## 3) data in memory as plot.data['step'], plot.data['fn'], plot.data['un'], etc. under the labels they were saved

