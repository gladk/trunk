#include"SimpleScene.hpp"
#include<yade/pkg-dem/Shop.hpp>
#include<yade/pkg-dem/BodyMacroParameters.hpp>
#include<yade/pkg-common/Bo1_Sphere_Aabb.hpp>
#include<yade/pkg-common/Bo1_Box_Aabb.hpp>
#include<yade/pkg-common/ForceResetter.hpp>
#include<yade/pkg-common/StateMetaEngine.hpp>
#include<yade/pkg-common/Aabb.hpp>
#include<yade/pkg-common/Box.hpp>
#include<yade/pkg-common/NewtonsForceLaw.hpp>
#include<yade/pkg-common/NewtonsMomentumLaw.hpp>
#include<yade/pkg-common/LeapFrogPositionIntegrator.hpp>
#include<yade/pkg-common/LeapFrogOrientationIntegrator.hpp>
#include<yade/pkg-dem/Ig2_Sphere_Sphere_ScGeom.hpp>
#include<yade/pkg-dem/Ig2_Box_Sphere_ScGeom.hpp>
#include<yade/pkg-common/StateMetaEngine.hpp>
#include<yade/pkg-common/InteractionGeometryDispatcher.hpp>
#include<yade/pkg-common/InteractionPhysicsDispatcher.hpp>
#include<yade/pkg-common/BoundDispatcher.hpp>
#include<yade/pkg-common/PhysicalActionDamper.hpp>
#include<yade/pkg-common/PhysicalActionApplier.hpp>
#include<yade/pkg-common/CundallNonViscousDamping.hpp>
#include<yade/pkg-common/CundallNonViscousDamping.hpp>
#include<yade/pkg-common/GravityEngines.hpp>
#include<yade/pkg-dem/ElasticContactLaw.hpp>
#include<yade/pkg-dem/ScGeom.hpp>
#include<yade/pkg-dem/Ip2_FrictMat_FrictMat_FrictPhys.hpp>
#include<yade/pkg-common/InsertionSortCollider.hpp>



YADE_PLUGIN((SimpleScene));
bool SimpleScene::generate(){
	message="";
	//@
	rootBody=Shop::rootBody();
	//@
	/* initializers */
		rootBody->initializers.clear();
		//@
		shared_ptr<BoundDispatcher> boundDispatcher	= shared_ptr<BoundDispatcher>(new BoundDispatcher);
			boundDispatcher->add(new Bo1_Sphere_Aabb);
			boundDispatcher->add(new Bo1_Box_Aabb);
			rootBody->initializers.push_back(boundDispatcher);
	//@
	/* engines */
		rootBody->engines.clear();
		//@
		rootBody->engines.push_back(shared_ptr<Engine>(new ForceResetter));
		//@
		// use boundDispatcher that we defined above
		rootBody->engines.push_back(boundDispatcher);
		//@
		shared_ptr<InsertionSortCollider> collider(new InsertionSortCollider);
			rootBody->engines.push_back(collider);
		//@
		shared_ptr<InteractionGeometryDispatcher> igeomDispatcher(new InteractionGeometryDispatcher);
			igeomDispatcher->add(new Ig2_Sphere_Sphere_ScGeom);
			igeomDispatcher->add(new Ig2_Box_Sphere_ScGeom);
			rootBody->engines.push_back(igeomDispatcher);
		//@
		shared_ptr<InteractionPhysicsDispatcher> iphysDispatcher(new InteractionPhysicsDispatcher);
			iphysDispatcher->add(new Ip2_FrictMat_FrictMat_FrictPhys);
			rootBody->engines.push_back(iphysDispatcher);
		//@
		shared_ptr<ElasticContactLaw> ecl(new ElasticContactLaw);
			rootBody->engines.push_back(ecl);
		//@
		shared_ptr<GravityEngine> ge(new GravityEngine);
			ge->gravity=Vector3r(0,0,-9.81);
			rootBody->engines.push_back(ge);
		//@
		shared_ptr<PhysicalActionDamper> dampingDispatcher(new PhysicalActionDamper);
			shared_ptr<CundallNonViscousForceDamping> forceDamper(new CundallNonViscousForceDamping);
			forceDamper->damping=0.2;
			dampingDispatcher->add(forceDamper);
			shared_ptr<CundallNonViscousMomentumDamping> momentumDamper(new CundallNonViscousMomentumDamping);
			momentumDamper->damping=0.2;
			dampingDispatcher->add(momentumDamper);
			rootBody->engines.push_back(dampingDispatcher);
		//@
		shared_ptr<PhysicalActionApplier> applyActionDispatcher(new PhysicalActionApplier);
			applyActionDispatcher->add(new NewtonsForceLaw);
			applyActionDispatcher->add(new NewtonsMomentumLaw);
			rootBody->engines.push_back(applyActionDispatcher);
		//@
		shared_ptr<StateMetaEngine> positionIntegrator(new StateMetaEngine);
			positionIntegrator->add(new LeapFrogPositionIntegrator);
			rootBody->engines.push_back(positionIntegrator);
		//@
		shared_ptr<StateMetaEngine> orientationIntegrator(new StateMetaEngine);
			orientationIntegrator->add(new LeapFrogOrientationIntegrator);
			rootBody->engines.push_back(orientationIntegrator);
	//@
	// set default values for Shop	
	Shop::setDefault("body_sdecGroupMask",1);
	Shop::setDefault("phys_density",2400.); Shop::setDefault("phys_young",30e9); Shop::setDefault("phys_poisson",.3);
	Shop::setDefault("aabb_randomColor",false);Shop::setDefault("shape_randomColor",false); Shop::setDefault("mold_randomColor",false);
	Shop::setDefault("aabb_color",Vector3r(0,1,0)); Shop::setDefault("shape_color",Vector3r(1,0,0)); Shop::setDefault("mold_color",Vector3r(1,0,0));

	//@
	shared_ptr<Body> box=Shop::box(Vector3r(0,0,0),Vector3r(.5,.5,.5));
	box->isDynamic=false;
	rootBody->bodies->insert(box);
	
	//@
	if(false){
		Vector3r extents(.5,.5,.5);
		shared_ptr<Body> b=shared_ptr<Body>(new Body(body_id_t(0),0));
		b->isDynamic=true;
		
		// phys
		shared_ptr<BodyMacroParameters> physics(new BodyMacroParameters);
		physics->mass=8*extents[0]*extents[1]*extents[2]*2400;
		physics->inertia=Vector3r(physics->mass*(4*extents[1]*extents[1]+4*extents[2]*extents[2])/12.,physics->mass*(4*extents[0]*extents[0]+4*extents[2]*extents[2])/12.,physics->mass*(4*extents[0]*extents[0]+4*extents[1]*extents[1])/12.);
		physics->se3=Se3r(Vector3r(0,0,0),Quaternionr::Identity());
		physics->young=30e9;
		physics->poisson=.3;
		b->physicalParameters=physics;

		// aabb
		shared_ptr<Aabb> aabb(new Aabb);
		aabb->diffuseColor=Vector3r(0,1,0);
		b->bound=aabb;

		// mold
		shared_ptr<Box> mold(new Box);
		mold->extents=extents;
		mold->diffuseColor=Vector3r(1,0,0);
		b->shape=mold;

		rootBody->bodies->insert(b);
	}

	//@
	Shop::setDefault("shape_color",Vector3r(0,1,0)); Shop::setDefault("mold_color",Vector3r(0,1,0));
	shared_ptr<Body> sphere(Shop::sphere(Vector3r(0,0,2),1));
	rootBody->bodies->insert(sphere);
	
	//@
	rootBody->dt=.2*Shop::PWaveTimeStep(rootBody);

	//@
	return true;
}

YADE_REQUIRE_FEATURE(PHYSPAR);

