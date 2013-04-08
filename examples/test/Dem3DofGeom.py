"Script showing shear interaction between facet/wall and sphere."
O.bodies.append([
	#sphere([0,0,0],1,dynamic=False,color=(0,1,0),wire=True),
	facet(([2,2,1],[-2,0,1],[2,-2,1]),fixed=True,color=(0,1,0),wire=False),
	#wall([0,0,1],axis=2,color=(0,1,0)),
	sphere([-1,0,2],1,fixed=False,color=(1,0,0),wire=True),
])
O.engines=[
	ForceResetter(),
	InsertionSortCollider([
		Bo1_Sphere_Aabb(),Bo1_Facet_Aabb(),Bo1_Wall_Aabb()
	]),
	IGeomDispatcher([
		Ig2_Sphere_Sphere_Dem3DofGeom(),
		Ig2_Facet_Sphere_Dem3DofGeom(),
		Ig2_Wall_Sphere_Dem3DofGeom()
	]),
	RotationEngine(rotationAxis=[0,1,0],angularVelocity=10,ids=[1]),
	TranslationEngine(translationAxis=[1,0,0],velocity=10,ids=[1]),
	NewtonIntegrator()#gravity=(0,0,-10))
]
O.miscParams=[
	Gl1_Dem3DofGeom_SphereSphere(normal=True,rolledPoints=True,unrolledPoints=True,shear=True,shearLabel=True),
	Gl1_Dem3DofGeom_FacetSphere(normal=False,rolledPoints=True,unrolledPoints=True,shear=True,shearLabel=True),
	Gl1_Dem3DofGeom_WallSphere(normal=False,rolledPoints=True,unrolledPoints=True,shear=True,shearLabel=True),
	Gl1_Sphere(wire=True)
]

try:
	from yade import qt
	renderer=qt.Renderer()
	renderer.wire=True
	renderer.intrGeom=True
	qt.Controller()
	qt.View()
except ImportError: pass

O.dt=1e-6
O.saveTmp('init')
O.run(1)
