# -*- encoding=utf-8 -*-

O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Box_Aabb()]),
	IGeomDispatcher([
		Ig2_Sphere_Sphere_ScGeom(),
		Ig2_Facet_Sphere_ScGeom(),
	]),
	IPhysDispatcher([Ip2_FrictMat_FrictMat_FrictPhys()]),
	ElasticContactLaw(),
	RotationEngine(ids=[1],rotationAxis=(1,0,0),angularVelocity=.01),
	RotationEngine(ids=[0],rotationAxis=(1,1,1),angularVelocity=-.02),
	NewtonIntegrator(damping=.2,gravity=(0,0,-9.81))
]

O.bodies.append(sphere([0,0,0],1,fixed=True,color=[1,0,0],wire=True))
O.bodies.append(sphere([0,sqrt(2),sqrt(2)],1,color=[0,1,0],wire=True))

O.dt=.01*PWaveTimeStep()
O.saveTmp()
#o.run(100000); o.wait(); print o.iter/o.realtime,'iterations/sec'
from yade import qt
qt.View()
renderer=qt.Renderer()
renderer.intrGeom=True
qt.Controller()
O.step(); O.step(); O.step()
O.run(20000)
