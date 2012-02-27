
#include<yade/core/Scene.hpp>
#include<yade/pkg-dem/Shop.hpp>

#include<yade/pkg-common/BoundDispatcher.hpp>
#include<yade/pkg-common/InteractionGeometryDispatcher.hpp>
#include<yade/pkg-common/ForceResetter.hpp>
#include<yade/pkg-common/InsertionSortCollider.hpp>
#include<yade/pkg-common/GravityEngines.hpp>
#include<yade/pkg-dem/NewtonIntegrator.hpp>

#include"TetraTestGen.hpp"

#include<yade/pkg-dem/Tetra.hpp>

YADE_PLUGIN((TetraTestGen));
bool TetraTestGen::generate()
{
	Shop::setDefault<int>("param_timeStepUpdateInterval",-1);

	rootBody=Shop::rootBody();

	shared_ptr<BoundDispatcher> boundDispatcher	= shared_ptr<BoundDispatcher>(new BoundDispatcher);
	boundDispatcher->add(new Bo1_Tetra_Aabb);
	rootBody->initializers.push_back(boundDispatcher);
	rootBody->engines.clear();
	rootBody->engines.push_back(shared_ptr<Engine>(new ForceResetter));
	rootBody->engines.push_back(boundDispatcher);
	rootBody->engines.push_back(shared_ptr<Engine>(new InsertionSortCollider));
	shared_ptr<InteractionGeometryDispatcher> interactionGeometryDispatcher(new InteractionGeometryDispatcher);
	interactionGeometryDispatcher->add(new Ig2_Tetra_Tetra_TTetraGeom);
	rootBody->engines.push_back(interactionGeometryDispatcher);
	// do not add any InteractionPhysicsDispatcher
	shared_ptr<TetraVolumetricLaw> constitutiveLaw(new TetraVolumetricLaw);
	rootBody->engines.push_back(constitutiveLaw);

	shared_ptr<GravityEngine> gravityCondition(new GravityEngine);
	gravityCondition->gravity=Vector3r(0,0,-9.81);
	rootBody->engines.push_back(gravityCondition);
	rootBody->engines.push_back(shared_ptr<Engine>(new NewtonIntegrator));

	rootBody->dt=1e-5;
	
	#if 0
	shared_ptr<Scene> oldRootBody=Omega::instance().getScene();
	Omega::instance().setScene(rootBody);
	#endif

	Vector3r v[4];

	v[0]=Vector3r(1,1,0);
	v[1]=Vector3r(0,-2,0);
	v[2]=Vector3r(-1,1,0);
	v[3]=Vector3r(0,0,-2);
	shared_ptr<Body> ground=Shop::tetra(v);
	ground->isDynamic=false;
	rootBody->bodies->insert(ground);


	for(size_t i=0; i<gridSize[0]; i++) for(size_t j=0; j<gridSize[1];j++) for(size_t k=0; k<gridSize[2];k++){
		Vector3r center=(1./std::max(gridSize[0],gridSize[1]))*Vector3r(i-gridSize[0]/2.,j-gridSize[1]/2.,.5+k);
		double size=.5/std::max(gridSize[0],gridSize[1]);
		#define __VERTEX(x,y,z) center+size*(Vector3r(x,y,z)+Vector3r(Mathr::IntervalRandom(-.5,.5),Mathr::IntervalRandom(-.5,.5),Mathr::IntervalRandom(-.5,.5)))
			v[0]=__VERTEX(1,1,0);
			v[1]=__VERTEX(0,-1,0);
			v[2]=__VERTEX(-1,0,0);
			v[3]=__VERTEX(0,0,2);
		#undef __VERTEX
		shared_ptr<Body> t=Shop::tetra(v);
		rootBody->bodies->insert(t);
	}

	// test inertia calculation
	#if 1
		v[0]=Vector3r(8.33220, -11.86875, 0.93355);
		v[1]=Vector3r(0.75523 ,5.00000, 16.37072);
		v[2]=Vector3r(52.61236, 5.00000, -5.38580);
		v[3]=Vector3r(2.00000, 5.00000, 3.00000);
		vector<Vector3r> vv;
		for(int i=0;i<4;i++)vv.push_back(v[i]);
		/* Vector3r cg=(v[0]+v[1]+v[2]+v[3])*.25;
		cerr<<"Centroid: "<<cg<<endl;
		v[0]-=cg; v[1]-=cg; v[2]-=cg; v[3]-=cg; */
		Matrix3r I=TetrahedronCentralInertiaTensor(vv);
		//cerr<<vv[0][0]<<" "<<vv[0][1]<<" "<<vv[0][2]<<endl;
		cerr<<I(0,0)<<endl<<I(1,1)<<endl<<I(2,2)<<endl<<-I(1,2)<<endl<<-I(0,1)<<endl<<-I(0,2)<<endl;
	#endif 

	#if 0
	Omega::instance().setScene(oldRootBody);
	#endif
	
	message="ATTN: this example is not working yet!";
	return true;
}

YADE_REQUIRE_FEATURE(PHYSPAR);

