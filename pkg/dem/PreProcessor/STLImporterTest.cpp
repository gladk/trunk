/*************************************************************************
*  Copyright (C) 2008 by Sergei Dorofeenko				 *
*  sega@users.berlios.de                                                 *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/


#include"STLImporterTest.hpp"
#include<yade/pkg-common/SpatialQuickSortCollider.hpp>
#include<yade/pkg-dem/STLImporter.hpp>
#include<yade/pkg-common/Facet.hpp>
#include<yade/core/Body.hpp>
#include<yade/core/Interaction.hpp>
#include<yade/core/Scene.hpp>
#include<yade/pkg-common/Aabb.hpp>
#include<yade/pkg-common/BoundDispatcher.hpp>
#include<yade/pkg-common/CundallNonViscousDamping.hpp>
#include<yade/pkg-common/CundallNonViscousDamping.hpp>
#include<yade/pkg-common/GravityEngines.hpp>
#include<yade/pkg-common/Sphere.hpp>
#include<yade/pkg-common/InteractionGeometryDispatcher.hpp>
#include<yade/pkg-common/InteractionPhysicsDispatcher.hpp>
#include<yade/pkg-common/PhysicalActionApplier.hpp>
#include<yade/pkg-common/ForceResetter.hpp>
#include<yade/pkg-common/PhysicalActionDamper.hpp>
#include<yade/pkg-common/StateMetaEngine.hpp>
#include<yade/pkg-dem/BodyMacroParameters.hpp>
#include<yade/pkg-dem/ElasticContactLaw.hpp>
#include<yade/pkg-dem/ElasticCriterionTimeStepper.hpp>
#include<yade/pkg-dem/MacroMicroElasticRelationships.hpp>
#include<yade/pkg-common/RotationEngine.hpp>

STLImporterTest::STLImporterTest() : FileGenerator()
{
	nbSpheres = Vector3r(10,10,10);
	minRadius = 0.07;
	maxRadius = 0.07;
	dampingForce = 0.3;
	dampingMomentum = 0.3;
	timeStepUpdateInterval = 300;
	sphereYoungModulus   = 15000000.0;
	spherePoissonRatio  = 0.2;
	sphereFrictionDeg   = 18.0;
	density = 2600;
	gravity = Vector3r(0,-9.81,0);
	disorder = Vector3r(0.002,0.002,0.002);
	wire=true;
	stlFileName = "hourglass.stl";
	angularVelocity = 0.5;
	rotationAxis = Vector3r(0,0,1);
	spheresHeight = 0;
}


STLImporterTest::~STLImporterTest()
{

}


void STLImporterTest::postProcessAttributes(bool)
{

}


bool STLImporterTest::generate()
{
	rootBody = shared_ptr<Scene>(new Scene);
	positionRootBody(rootBody);

	rootBody->dt = 0.001; //default time step

////////////////////////////////////
///////// Container
	
		
////////////////////////////////////
///////// walls

	cerr << "Import walls geometry from file: " << stlFileName << endl;
	STLImporter imp;
	if (!imp.open(stlFileName.c_str()))
	{
	    cerr << "ERROR: Bad file: " << stlFileName << endl;
	    message="Input file not found, you can copy it from examples/ directory or make one using blender 3D modelling";
	    return false;
	}
	imp.wire=wire;
	// create bodies
	for(int i=0,e=imp.number_of_facets;i<e;++i)
	{
	    shared_ptr<Body> b(new Body(body_id_t(0),1));
    
	    b->isDynamic		= false;
	    
	    // physical parameters
	    shared_ptr<BodyMacroParameters> physics(new BodyMacroParameters);
	    physics->angularVelocity	= Vector3r(0,0,0);
	    physics->velocity		= Vector3r(0,0,0);
	    physics->mass			= 0;
	    physics->inertia		= Vector3r(0,0,0);
	    physics->young			= sphereYoungModulus;
	    physics->poisson		= spherePoissonRatio;
	    physics->frictionAngle		= sphereFrictionDeg * Mathr::PI/180.0;
	    b->physicalParameters	= physics;

	    // bounding box 
		shared_ptr<Aabb> aabb(new Aabb);
		aabb->diffuseColor		= Vector3r(0,1,0);
		b->bound	= aabb;
	    
	    rootBody->bodies->insert(b);
	}
	// import bodies (create geometry)
	imp.import(rootBody->bodies);

///////// spheres
	float all = nbSpheres[0]*nbSpheres[1]*nbSpheres[2];
	float current = 0.0;

	for(int i=0;i<nbSpheres[0];i++)
	{
		if(shouldTerminate()) return false;

		for(int j=0;j<nbSpheres[1];j++)
			for(int k=0;k<nbSpheres[2];k++)
			{
				shared_ptr<Body> sphere;
				createSphere(sphere,i,j,k);
				rootBody->bodies->insert(sphere);

				setProgress(current++/all);
			}
	}
	

///////// engines
	createActors(rootBody);

	return true;
}

void STLImporterTest::createSphere(shared_ptr<Body>& body, int i, int j, int k)
{
	body = shared_ptr<Body>(new Body(body_id_t(0),1));
	shared_ptr<BodyMacroParameters> physics(new BodyMacroParameters);
	shared_ptr<Aabb> aabb(new Aabb);
	shared_ptr<Sphere> iSphere(new Sphere);
	
	Vector3r position		= Vector3r(i,j+spheresHeight,k)*(2*maxRadius*1.1) 
					  - Vector3r( nbSpheres[0]/2*(2*maxRadius*1.1) , 0 , nbSpheres[2]/2*(2*maxRadius*1.1) )
					  + Vector3r( 	 Mathr::SymmetricRandom()*disorder[0]
					  		,Mathr::SymmetricRandom()*disorder[1]
							,Mathr::SymmetricRandom()*disorder[2])*maxRadius;
	
	Real radius 			= (Mathr::IntervalRandom(minRadius,maxRadius));
	
	body->isDynamic			= true;
	
	physics->angularVelocity	= Vector3r(0,0,0);
	physics->velocity		= Vector3r(0,0,0);
	physics->mass			= 4.0/3.0*Mathr::PI*radius*radius*radius*density;
	physics->inertia		= Vector3r(2.0/5.0*physics->mass*radius*radius,2.0/5.0*physics->mass*radius*radius,2.0/5.0*physics->mass*radius*radius); //
	physics->se3			= Se3r(position,Quaternionr::Identity());
	physics->young			= sphereYoungModulus;
	physics->poisson		= spherePoissonRatio;
	physics->frictionAngle		= sphereFrictionDeg * Mathr::PI/180.0;

	aabb->diffuseColor		= Vector3r(0,1,0);

	iSphere->radius			= radius;
	iSphere->diffuseColor		= Vector3r(0.8,0.3,0.3);

	body->shape	= iSphere;
	body->bound		= aabb;
	body->physicalParameters	= physics;
}

void STLImporterTest::createActors(shared_ptr<Scene>& rootBody)
{
	
	shared_ptr<InteractionGeometryDispatcher> interactionGeometryDispatcher(new InteractionGeometryDispatcher);
	interactionGeometryDispatcher->add("Ig2_Sphere_Sphere_ScGeom");
	interactionGeometryDispatcher->add("Ig2_Facet_Sphere_ScGeom");

	shared_ptr<InteractionPhysicsDispatcher> interactionPhysicsDispatcher(new InteractionPhysicsDispatcher);
	interactionPhysicsDispatcher->add("MacroMicroElasticRelationships");
		
	shared_ptr<BoundDispatcher> boundDispatcher	= shared_ptr<BoundDispatcher>(new BoundDispatcher);
	boundDispatcher->add("Bo1_Sphere_Aabb");
	boundDispatcher->add("Bo1_Facet_Aabb");
	
	shared_ptr<GravityEngine> gravityCondition(new GravityEngine);
	gravityCondition->gravity = gravity;
	
	shared_ptr<CundallNonViscousForceDamping> actionForceDamping(new CundallNonViscousForceDamping);
	actionForceDamping->damping = dampingForce;
	shared_ptr<CundallNonViscousMomentumDamping> actionMomentumDamping(new CundallNonViscousMomentumDamping);
	actionMomentumDamping->damping = dampingMomentum;
	shared_ptr<PhysicalActionDamper> actionDampingDispatcher(new PhysicalActionDamper);
	actionDampingDispatcher->add(actionForceDamping);
	actionDampingDispatcher->add(actionMomentumDamping);
	
	shared_ptr<PhysicalActionApplier> applyActionDispatcher(new PhysicalActionApplier);
	applyActionDispatcher->add("NewtonsForceLaw");
	applyActionDispatcher->add("NewtonsMomentumLaw");
	
	shared_ptr<StateMetaEngine> positionIntegrator(new StateMetaEngine);
	positionIntegrator->add("LeapFrogPositionIntegrator");
	shared_ptr<StateMetaEngine> orientationIntegrator(new StateMetaEngine);
	orientationIntegrator->add("LeapFrogOrientationIntegrator");

	shared_ptr<RotationEngine> kinematic = shared_ptr<RotationEngine>(new RotationEngine);
 	kinematic->angularVelocity  = angularVelocity;
	rotationAxis.normalize();
 	kinematic->rotationAxis  = rotationAxis;
 	kinematic->rotateAroundZero = true;
	
	
	shared_ptr<Shape> facet(new Facet);
	for(BodyContainer::iterator bi = rootBody->bodies->begin(), biEnd=rootBody->bodies->end(); bi!=biEnd; ++bi)
	    if ( (*bi)->shape->getClassIndex() == facet->getClassIndex() )
		kinematic->subscribedBodies.push_back((*bi)->getId());

	shared_ptr<ElasticCriterionTimeStepper> sdecTimeStepper(new ElasticCriterionTimeStepper);
	sdecTimeStepper->timeStepUpdateInterval = timeStepUpdateInterval;

	rootBody->engines.clear();
	rootBody->engines.push_back(shared_ptr<Engine>(new ForceResetter));
	rootBody->engines.push_back(sdecTimeStepper);
	rootBody->engines.push_back(boundDispatcher);	
	rootBody->engines.push_back(shared_ptr<Engine>(new SpatialQuickSortCollider));
	rootBody->engines.push_back(interactionGeometryDispatcher);
	rootBody->engines.push_back(interactionPhysicsDispatcher);
	rootBody->engines.push_back(shared_ptr<Engine>(new ElasticContactLaw));
	rootBody->engines.push_back(gravityCondition);
	rootBody->engines.push_back(actionDampingDispatcher);
	rootBody->engines.push_back(applyActionDispatcher);
	rootBody->engines.push_back(positionIntegrator);
	rootBody->engines.push_back(orientationIntegrator);
	rootBody->engines.push_back(kinematic);
 	
	rootBody->initializers.clear();
	rootBody->initializers.push_back(boundDispatcher);
}


void STLImporterTest::positionRootBody(shared_ptr<Scene>& rootBody) 
{
	rootBody->isDynamic		= false;
	
	shared_ptr<ParticleParameters> physics(new ParticleParameters); // FIXME : fix indexable class PhysicalParameters
	physics->se3				= Se3r(Vector3r(0,0,0),Quaternionr::Identity());
	physics->mass				= 0;
	physics->velocity			= Vector3r(0,0,0);
	physics->acceleration			= Vector3r::Zero();
		
	shared_ptr<Aabb> aabb(new Aabb);
	aabb->diffuseColor			= Vector3r(0,0,1);
	
	rootBody->bound		= YADE_PTR_CAST<Bound>(aabb);
	rootBody->physicalParameters 		= physics;
}


YADE_PLUGIN((STLImporterTest));

YADE_REQUIRE_FEATURE(PHYSPAR);

