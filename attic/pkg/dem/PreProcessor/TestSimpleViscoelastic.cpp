/*************************************************************************
*  Copyright (C) 2006 by Sergei Dorofeenko				 *
*  sega@users.berlios.de                                                 *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#include "TestSimpleViscoelastic.hpp"
#include <yade/pkg-dem/Shop.hpp>
#include <yade/core/Body.hpp>
#include <yade/core/Scene.hpp>
#include<yade/core/Body.hpp>
#include<yade/core/Interaction.hpp>
#include<yade/core/Scene.hpp>
#include<yade/pkg-common/Aabb.hpp>
#include<yade/pkg-common/BoundDispatcher.hpp>
#include<yade/pkg-common/Sphere.hpp>
#include<yade/pkg-common/Box.hpp>
#include<yade/pkg-common/InteractionGeometryDispatcher.hpp>
#include<yade/pkg-common/InteractionPhysicsDispatcher.hpp>
#include<yade/pkg-common/InsertionSortCollider.hpp>
#include<yade/pkg-common/PhysicalActionApplier.hpp>
#include<yade/pkg-common/ForceResetter.hpp>
#include<yade/pkg-common/StateMetaEngine.hpp>
#include<yade/pkg-common/SphereModel.hpp>
#include<yade/pkg-common/BoxModel.hpp>
#include<yade/pkg-common/LawDispatcher.hpp>
#include<yade/pkg-dem/RigidBodyRecorder.hpp>
#include<yade/pkg-dem/SimpleViscoelasticSpheresInteractionRecorder.hpp>
#include<yade/pkg-dem/SimpleViscoelasticBodyParameters.hpp>
#include<yade/pkg-dem/ef2_Spheres_Viscoelastic_SimpleViscoelasticContactLaw.hpp>
#include<yade/pkg-dem/SimpleViscoelasticRelationships.hpp>
#include<yade/pkg-common/GravityEngines.hpp>

YADE_REQUIRE_FEATURE(geometricalmodel)

TestSimpleViscoelastic::TestSimpleViscoelastic() : FileGenerator()
{
    tc = 0.001;
    en = 0.3;
    es = 0.3;
    frictionAngle = 0.52; // 30 grad
    R  = 1;
    density = 2600;
    nbSpheres = 10;
    h  = 1;
    velocity = Vector3r(0,-1,0);
    angular_vel = Vector3r(0,0,0);
    rotationBlocked=false;
    groundSize = Vector3r(10,1,10);
    gravity = Vector3r(0,-9.81,0);
    outputBase="test";
    interval=50;
}

TestSimpleViscoelastic::~TestSimpleViscoelastic()
{}

void TestSimpleViscoelastic::postProcessAttributes(bool)
{}

bool TestSimpleViscoelastic::generate()
{
    rootBody = shared_ptr<Scene>(new Scene);
    createActors(rootBody);
    positionRootBody(rootBody);
    
////////////////////////////////////
///////// Container
	
////////////////////////////////////
///////// ground
    
    shared_ptr<Body> ground;
    createBox(ground, Vector3r(0,0,0),groundSize);
    rootBody->bodies->insert(ground);

////////////////////////////////////
///////// Spheres

    float all = nbSpheres;
    float current = 0.0;
    for(int i=0; i<nbSpheres; ++i)
    {
	if (shouldTerminate()) return false;
	shared_ptr<Body> sphere;
	createSphere(sphere,i);
	rootBody->bodies->insert(sphere);

	setProgress(current++/all);
    }
    return true;
}

void TestSimpleViscoelastic::createBox(shared_ptr<Body>& body, Vector3r position, Vector3r extents)
{
	body = shared_ptr<Body>(new Body(body_id_t(0),1));
	shared_ptr<SimpleViscoelasticBodyParameters> physics(new SimpleViscoelasticBodyParameters);
	shared_ptr<Aabb> aabb(new Aabb);
	shared_ptr<BoxModel> gBox(new BoxModel);
	shared_ptr<Box> iBox(new Box);
	
	
	Quaternionr q;
	q.FromAxisAngle( Vector3r(0,0,1),0);

	body->isDynamic			= false;
	
	physics->angularVelocity	= Vector3r(0,0,0);
	physics->velocity		= Vector3r(0,0,0);
	physics->mass			= extents[0]*extents[1]*extents[2]*density*8; 
	physics->inertia		= Vector3r(physics->mass*(extents[1]*extents[1]+extents[2]*extents[2])/3
						, physics->mass*(extents[0]*extents[0]+extents[2]*extents[2])/3
						, physics->mass*(extents[1]*extents[1]+extents[0]*extents[0])/3);
	physics->se3			= Se3r(position,q);
    Shop::getViscoelasticFromSpheresInteraction(physics->mass,tc,en,es,physics);
	physics->frictionAngle			= frictionAngle;

	aabb->diffuseColor		= Vector3r(1,0,0);

	gBox->extents			= extents;
	gBox->diffuseColor		= Vector3r(1,1,1);
	gBox->wire			= false;
	gBox->shadowCaster		= true;
	
	iBox->extents			= extents;
	iBox->diffuseColor		= Vector3r(1,1,1);

	body->bound		= aabb;
	body->shape	= iBox;
	body->geometricalModel		= gBox;
	body->physicalParameters	= physics;
}

void TestSimpleViscoelastic::createActors(shared_ptr<Scene>& rootBody)
{
    
    shared_ptr<InteractionGeometryDispatcher> interactionGeometryDispatcher(new InteractionGeometryDispatcher);
    interactionGeometryDispatcher->add("Ig2_Sphere_Sphere_ScGeom");
    interactionGeometryDispatcher->add("Ig2_Box_Sphere_ScGeom");

    shared_ptr<InteractionPhysicsDispatcher> interactionPhysicsDispatcher(new InteractionPhysicsDispatcher);
    interactionPhysicsDispatcher->add("SimpleViscoelasticRelationships");
	    
    shared_ptr<BoundDispatcher> boundDispatcher	= shared_ptr<BoundDispatcher>(new BoundDispatcher);
    boundDispatcher->add("Bo1_Sphere_Aabb");
    boundDispatcher->add("Bo1_Box_Aabb");
    
    shared_ptr<GravityEngine> gravityCondition(new GravityEngine);
    gravityCondition->gravity = gravity;
    
    shared_ptr<PhysicalActionApplier> applyActionDispatcher(new PhysicalActionApplier);
    applyActionDispatcher->add("NewtonsForceLaw");
    applyActionDispatcher->add("NewtonsMomentumLaw");
    
    shared_ptr<StateMetaEngine> positionIntegrator(new StateMetaEngine);
    positionIntegrator->add("LeapFrogPositionIntegrator");
    shared_ptr<StateMetaEngine> orientationIntegrator(new StateMetaEngine);
    orientationIntegrator->add("LeapFrogOrientationIntegrator");
    
	shared_ptr<LawDispatcher> constitutiveLaw(new LawDispatcher);
	constitutiveLaw->add("ef2_Spheres_Viscoelastic_SimpleViscoelasticContactLaw");

    rootBody->engines.clear();
    rootBody->engines.push_back(shared_ptr<Engine>(new ForceResetter));
    rootBody->engines.push_back(boundDispatcher);	
    rootBody->engines.push_back(shared_ptr<Engine>(new InsertionSortCollider));
    rootBody->engines.push_back(interactionGeometryDispatcher);
    rootBody->engines.push_back(interactionPhysicsDispatcher);
    rootBody->engines.push_back(constitutiveLaw);
    rootBody->engines.push_back(gravityCondition);
    rootBody->engines.push_back(applyActionDispatcher);
    rootBody->engines.push_back(positionIntegrator);
    if(!rotationBlocked) rootBody->engines.push_back(orientationIntegrator);

    shared_ptr<RigidBodyRecorder> rigidBodyRecorder(new RigidBodyRecorder);
    rigidBodyRecorder->outputBase = outputBase+string("-body");
    rigidBodyRecorder->interval = interval;
    rootBody->engines.push_back(rigidBodyRecorder);

    shared_ptr<SimpleViscoelasticSpheresInteractionRecorder> interactionRecorder(new SimpleViscoelasticSpheresInteractionRecorder);
    interactionRecorder->outputBase = outputBase+string("-interaction");
    interactionRecorder->interval = interval;
    rootBody->engines.push_back(interactionRecorder);

    rootBody->initializers.clear();
    rootBody->initializers.push_back(boundDispatcher);
}

void TestSimpleViscoelastic::createSphere(shared_ptr<Body>& body, int i)
{
	body = shared_ptr<Body>(new Body(body_id_t(0),1));
	shared_ptr<SimpleViscoelasticBodyParameters> physics(new SimpleViscoelasticBodyParameters);
	shared_ptr<Aabb> aabb(new Aabb);
	shared_ptr<SphereModel> gSphere(new SphereModel);
	shared_ptr<Sphere> iSphere(new Sphere);
	
	Quaternionr q;
	q.FromAxisAngle( Vector3r(0,0,1),0);
	
	Vector3r position = Vector3r(0,i*(2.0*R+h)+2*groundSize[1]+h,0);
	
//	Real radius 			= (Mathr::IntervalRandom(minRadius,maxRadius));
	Real radius = R;
	
	body->isDynamic	= true;
	
	physics->angularVelocity	= angular_vel;
	physics->velocity		= velocity;
	physics->mass			= 4.0/3.0*Mathr::PI*radius*radius*radius*density;
	physics->inertia		= Vector3r(2.0/5.0*physics->mass*radius*radius,2.0/5.0*physics->mass*radius*radius,2.0/5.0*physics->mass*radius*radius); //
	physics->se3			= Se3r(position,q);
    Shop::getViscoelasticFromSpheresInteraction(physics->mass,tc,en,es,physics);
    physics->frictionAngle		= frictionAngle;

	aabb->diffuseColor		= Vector3r(0,1,0);

	gSphere->radius			= radius;
	gSphere->diffuseColor		= Vector3r(Mathr::UnitRandom(),Mathr::UnitRandom(),Mathr::UnitRandom());
	gSphere->wire			= false;
	gSphere->shadowCaster		= true;
	
	iSphere->radius			= radius;
	iSphere->diffuseColor		= Vector3r(0.8,0.3,0.3);

	body->shape	= iSphere;
	body->geometricalModel		= gSphere;
	body->bound		= aabb;
	body->physicalParameters	= physics;
}
    

void TestSimpleViscoelastic::positionRootBody(shared_ptr<Scene>& rootBody) 
{
	rootBody->isDynamic		= false;
	
	Quaternionr q;
	q.FromAxisAngle( Vector3r(0,0,1),0);

	shared_ptr<ParticleParameters> physics(new ParticleParameters); 
	physics->se3				= Se3r(Vector3r(0,0,0),q);
	physics->mass				= 0;
	physics->velocity			= Vector3r(0,0,0);
	physics->acceleration			= Vector3r::Zero();
		
	shared_ptr<Aabb> aabb(new Aabb);
	aabb->diffuseColor			= Vector3r(0,0,1);
	
	rootBody->bound		= YADE_PTR_CAST<Bound>(aabb);
	rootBody->physicalParameters 		= physics;
}

    
    

YADE_REQUIRE_FEATURE(PHYSPAR);

