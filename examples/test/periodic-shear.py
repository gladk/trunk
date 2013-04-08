O.periodic=True
O.cell.setBox(.55,.55,.55)
O.bodies.append(facet([[.4,.0001,.3],[.2,.0001,.3],[.3,.2,.2]]))
O.bodies.append(sphere([.3,.1,.4],.05,fixed=False))
O.bodies.append(sphere([.200001,.2000001,.4],.05,fixed=True))
O.bodies.append(sphere([.3,0,0],.1,fixed=True))
O.engines=[
	ForceResetter(),
	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Facet_Aabb()],label='isc'),
	InteractionLoop(
		[Ig2_Sphere_Sphere_ScGeom(),Ig2_Facet_Sphere_ScGeom()],
		[Ip2_FrictMat_FrictMat_FrictPhys()],
		[Law2_ScGeom_FrictPhys_CundallStrack()]
	),
	NewtonIntegrator(gravity=(0,0,-10)),
	#PyRunner(command='doCellFlip()',realPeriod=5)
]

def doCellFlip():
	flip=1 if O.cell.trsf[1,2]<0 else -1;
	flipCell(Matrix3(0,0,0, 0,0,flip, 0,0,0))

#g=0.
#while False:
#	O.cellShear=Vector3(.2*sin(g),.2*cos(pi*g),.2*sin(2*g)+.2*cos(3*g))
#	time.sleep(0.001)
#	g+=1e-3
O.cell.trsf=Matrix3(1,0,0, 0,1,.5, 0,0,1)
O.dt=2e-2*PWaveTimeStep()
O.step()
O.saveTmp()
rrr=yade.qt.Renderer()
rrr.intrAllWire,rrr.bound=True,True
#isc.watch1,isc.watch2=0,-1

#import yade.qt,time
#v=yade.qt.View()
#v.axes=True
#v.grid=(True,True,True)


