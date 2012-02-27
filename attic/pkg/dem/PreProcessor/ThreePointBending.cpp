/*************************************************************************
*  Copyright (C) 2007 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

///
// FIXME : add CohesiveContactLaw engine to make it work
///


#include "ThreePointBending.hpp"

#include<yade/pkg-dem/ElasticContactLaw.hpp>
#include<yade/pkg-dem/ElasticCohesiveLaw.hpp>
#include<yade/pkg-dem/MacroMicroElasticRelationships.hpp>
#include<yade/pkg-dem/BodyMacroParameters.hpp>
#include<yade/pkg-dem/ScGeom.hpp>
#include<yade/pkg-dem/SDECLinkPhysics.hpp>
#include<yade/pkg-dem/ElasticCriterionTimeStepper.hpp>


#include<yade/pkg-common/BoxModel.hpp>
#include<yade/pkg-common/Aabb.hpp>
#include<yade/pkg-common/SphereModel.hpp>
#include<yade/core/Scene.hpp>
#include<yade/pkg-common/InsertionSortCollider.hpp>
#include<yade/lib-serialization/IOFormatManager.hpp>
#include<yade/core/Interaction.hpp>
#include<yade/pkg-common/BoundDispatcher.hpp>

#include<yade/pkg-common/PhysicalActionDamper.hpp>
#include<yade/pkg-common/PhysicalActionApplier.hpp>
#include<yade/pkg-common/CundallNonViscousDamping.hpp>
#include<yade/pkg-common/CundallNonViscousDamping.hpp>

#include<yade/pkg-common/InteractionGeometryDispatcher.hpp>
#include<yade/pkg-common/InteractionPhysicsDispatcher.hpp>
#include<yade/core/Body.hpp>
#include<yade/pkg-common/Box.hpp>
#include<yade/pkg-common/Sphere.hpp>

#include<yade/pkg-common/ForceResetter.hpp>
#include<yade/pkg-common/StateMetaEngine.hpp>

#include<yade/pkg-common/TranslationEngine.hpp>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>

#include<yade/pkg-dem/Shop.hpp>

YADE_REQUIRE_FEATURE(geometricalmodel);

ThreePointBending::ThreePointBending () : FileGenerator()
{
	yadeFileWithSpheres	="";

	pistonVelocity		= 1.;
	dampingForce		= 0.2;
	dampingMomentum		= 0.2;
	
	timeStepUpdateInterval	= 300;
	momentRotationLaw       = true;
	
	//sphereYoungModulus	= 10000000;
	sphereYoungModulus	= 30000000;
	spherePoissonRatio	= 0.2;
	sphereFrictionDeg       = 18.0;
	
	linkKn			= 50000000;
	linkKs			= 5000000;
	linkMaxNormalForce	= 5500000;
	linkMaxShearForce       = 550000;
}


ThreePointBending::~ThreePointBending ()
{

}


void ThreePointBending::postProcessAttributes(bool)
{
}



bool ThreePointBending::generate()
{
	rootBody = shared_ptr<Scene>(new Scene);
	createActors(rootBody);
	positionRootBody(rootBody);

////////////////////////////////////
	
	Vector3r min(10000,10000,10000),max(-10000,-10000,-10000);

	// load simulation file, extract spheres and use those
	if (yadeFileWithSpheres.size()!=0){
		shared_ptr<Scene> metaBodyWithSpheres;
		IOFormatManager::loadFromFile("XMLFormatManager",yadeFileWithSpheres,"scene",metaBodyWithSpheres);
		assert(metaBodyWithSpheres->getClassName()=="Scene");

		FOREACH(shared_ptr<Body> b, *metaBodyWithSpheres->bodies){
			if(b->geometricalModel->getClassName()!="SphereModel") continue;
			min = componentMinVector(min,b->physicalParameters->se3.position - static_cast<SphereModel*>(b->geometricalModel.get())->radius * Vector3r(1,1,1));
			max = componentMaxVector(max,b->physicalParameters->se3.position + static_cast<SphereModel*>(b->geometricalModel.get())->radius * Vector3r(1,1,1));
			BodyMacroParameters* bm = dynamic_cast<BodyMacroParameters*>(b->physicalParameters.get());
			if(!bm) {message="Error: spheres don't use BodyMacroParameters for physical parameters"; return false;}
			bm->young		= sphereYoungModulus;
			bm->poisson		= spherePoissonRatio;
			bm->frictionAngle	= sphereFrictionDeg * Mathr::PI/180.0;
			bm->mass		/= 1000.0; // ???!

			rootBody->bodies->insert(b);
		}
	}
	// generate specimen
	else{
		Real radius=1; int xxSup=2, yySup=10, zzSup=2;
		min=Vector3r(-radius,-radius,-radius); max=Vector3r(radius*(2*xxSup-1),radius*(2*yySup-1),radius*(2*zzSup-1));
		for(int xx=0; xx<xxSup; xx++){
			for(int yy=0; yy<yySup; yy++){
				for(int zz=0; zz<zzSup; zz++){
					shared_ptr<Body> b=Shop::sphere(Vector3r(xx*2*radius,yy*2*radius,zz*2*radius),radius);
					rootBody->bodies->insert(b);
				}
			}
		}
	}
	
/////////////////////////////////////

	rootBody->interactions->clear();
	
	shared_ptr<Body> bodyA;

	BodyContainer::iterator bi    = rootBody->bodies->begin();
	BodyContainer::iterator biEnd = rootBody->bodies->end();
	BodyContainer::iterator bi2;

	Real interactionFactor=1.1;
	for( ; bi!=biEnd ; ++bi )
	{
		bodyA =*bi;
		bi2=bi;
		++bi2;
		for( ; bi2!=biEnd ; ++bi2 )
		{
			shared_ptr<Body> bodyB = *bi2;

			shared_ptr<BodyMacroParameters> a = YADE_PTR_CAST<BodyMacroParameters>(bodyA->physicalParameters);
			shared_ptr<BodyMacroParameters> b = YADE_PTR_CAST<BodyMacroParameters>(bodyB->physicalParameters);
			shared_ptr<Sphere>	as = YADE_PTR_CAST<Sphere>(bodyA->shape);
			shared_ptr<Sphere>	bs = YADE_PTR_CAST<Sphere>(bodyB->shape);

			if ((a->se3.position-b->se3.position).Length() < interactionFactor*(as->radius+bs->radius))  
			{
				shared_ptr<Interaction> 		link(new Interaction( bodyA->getId() , bodyB->getId() ));
				shared_ptr<ScGeom>		geometry(new ScGeom);
				shared_ptr<SDECLinkPhysics>	physics(new SDECLinkPhysics);
				
				geometry->radius1			= as->radius - fabs(as->radius - bs->radius)*0.5;
				geometry->radius2			= bs->radius - fabs(as->radius - bs->radius)*0.5;

				physics->initialKn			= linkKn; // FIXME - BIG problem here.
				physics->initialKs			= linkKs;
				physics->heta				= 1;
				physics->initialEquilibriumDistance	= (a->se3.position - b->se3.position).Length();
				physics->knMax				= linkMaxNormalForce;
				physics->ksMax				= linkMaxShearForce;

				link->interactionGeometry 		= geometry;
				link->interactionPhysics 		= physics;
				rootBody->interactions->insert(link);
			}
		}
	}

	// specimen supports
	shared_ptr<Body> piston;
	shared_ptr<Body> supportBox1;
	shared_ptr<Body> supportBox2;

	Real dimX=max[0]-min[0], dimY=max[1]-min[1], dimZ=max[2]-min[2];
	Vector3r center=.5*(max+min);

	createBox(piston,
		Vector3r(center[0],center[1],max[2]+.5*dimZ),
		Vector3r(.75*dimX,.05*dimY,.5*dimZ));
	createBox(supportBox1,
		Vector3r(center[0],min[1],min[2]-.5*dimZ),
		Vector3r(.5*dimX,.1*dimY,.5*dimZ));
	createBox(supportBox2,
		Vector3r(center[0],max[1],min[2]-.5*dimZ),
		Vector3r(.5*dimX,.1*dimY,.5*dimZ));
			
	rootBody->bodies->insert(piston);
	rootBody->bodies->insert(supportBox1);
	rootBody->bodies->insert(supportBox2);

	FOREACH(const shared_ptr<Engine>& e, rootBody->engines){
		if(e->getClassName()!="TranslationEngine") continue;
		shared_ptr<TranslationEngine> te=YADE_PTR_CAST<TranslationEngine>(e);
		te->subscribedBodies.clear(); te->subscribedBodies.push_back(piston->getId());
	}

	
	message="total number of permament links created: " 
		+ lexical_cast<string>(rootBody->interactions->size()) 
		+ "\nWARNING: link bonds are nearly working, but the formulas are waiting for total rewrite!"+
		+"\nWARNING: The results are meaningless, since ElasticCohesiveLaw works only with (unused) SDECLinkGeometry.";
	return true;
}


void ThreePointBending::createBox(shared_ptr<Body>& body, Vector3r position, Vector3r extents)
{
	body = shared_ptr<Body>(new Body(body_id_t(0),55));
	shared_ptr<BodyMacroParameters> physics(new BodyMacroParameters);
	shared_ptr<Aabb> aabb(new Aabb);
	shared_ptr<BoxModel> gBox(new BoxModel);
	shared_ptr<Box> iBox(new Box);
	
	body->isDynamic			= false;
	
	physics->angularVelocity	= Vector3r(0,0,0);
	physics->velocity		= Vector3r(0,0,0);
	physics->mass			= extents[0]*extents[1]*extents[2]*2; 
	physics->inertia		= Vector3r(
							  physics->mass*(extents[1]*extents[1]+extents[2]*extents[2])/3
							, physics->mass*(extents[0]*extents[0]+extents[2]*extents[2])/3
							, physics->mass*(extents[1]*extents[1]+extents[0]*extents[0])/3
						);
	physics->se3			= Se3r(position,Quaternionr::Identity());
	physics->young			= sphereYoungModulus;
	physics->poisson		= spherePoissonRatio;
	physics->frictionAngle		= sphereFrictionDeg * Mathr::PI/180.0;

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


void ThreePointBending::createActors(shared_ptr<Scene>& rootBody)
{
	
	shared_ptr<InteractionGeometryDispatcher> interactionGeometryDispatcher(new InteractionGeometryDispatcher);
	interactionGeometryDispatcher->add("Ig2_Sphere_Sphere_ScGeom");
	interactionGeometryDispatcher->add("Ig2_Box_Sphere_ScGeom");

	shared_ptr<InteractionPhysicsDispatcher> interactionPhysicsDispatcher(new InteractionPhysicsDispatcher);
	interactionPhysicsDispatcher->add("MacroMicroElasticRelationships");
		
	shared_ptr<BoundDispatcher> boundDispatcher	= shared_ptr<BoundDispatcher>(new BoundDispatcher);
	boundDispatcher->add("Bo1_Sphere_Aabb");
	boundDispatcher->add("Bo1_Box_Aabb");
		
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
 	
	shared_ptr<ElasticCriterionTimeStepper> sdecTimeStepper(new ElasticCriterionTimeStepper);
	sdecTimeStepper->sdecGroupMask = 55;
	sdecTimeStepper->timeStepUpdateInterval = timeStepUpdateInterval;

// moving wall
	shared_ptr<TranslationEngine> kinematic = shared_ptr<TranslationEngine>(new TranslationEngine);
	kinematic->velocity  = pistonVelocity;
	kinematic->translationAxis  = Vector3r(0,0,-1);
	// subscribed bodies will be initialized later

	shared_ptr<ElasticContactLaw> constitutiveLaw(new ElasticContactLaw);
	constitutiveLaw->sdecGroupMask = 55;
	constitutiveLaw->momentRotationLaw = momentRotationLaw;
// FIXME FIXME FIXME ....	
	shared_ptr<ElasticCohesiveLaw> constitutiveLaw2(new ElasticCohesiveLaw);
	constitutiveLaw2->sdecGroupMask = 55;
	constitutiveLaw2->momentRotationLaw = momentRotationLaw;
	
	rootBody->engines.clear();
	rootBody->engines.push_back(sdecTimeStepper);
	rootBody->engines.push_back(shared_ptr<Engine>(new ForceResetter));
	rootBody->engines.push_back(boundDispatcher);
	rootBody->engines.push_back(shared_ptr<Engine>(new InsertionSortCollider));
	rootBody->engines.push_back(interactionGeometryDispatcher);
	rootBody->engines.push_back(interactionPhysicsDispatcher);
	rootBody->engines.push_back(constitutiveLaw);
	rootBody->engines.push_back(constitutiveLaw2);
	rootBody->engines.push_back(actionDampingDispatcher);
	rootBody->engines.push_back(applyActionDispatcher);
	rootBody->engines.push_back(positionIntegrator);
	rootBody->engines.push_back(orientationIntegrator);
	rootBody->engines.push_back(kinematic);

	rootBody->initializers.clear();
	rootBody->initializers.push_back(boundDispatcher);
}
	

void ThreePointBending::positionRootBody(shared_ptr<Scene>& rootBody)
{
	rootBody->isDynamic		= false;
	
	Quaternionr q;
	q.FromAxisAngle( Vector3r(0,0,1),0);
	
	shared_ptr<ParticleParameters> physics(new ParticleParameters); // FIXME : fix indexable class PhysicalParameters
	physics->se3			= Se3r(Vector3r(0,0,0),q);
	physics->mass			= 0;
	physics->velocity		= Vector3r::Zero();
	physics->acceleration		= Vector3r::Zero();
	
	shared_ptr<Aabb> aabb(new Aabb);
	aabb->diffuseColor		= Vector3r(0,0,1);
	
	rootBody->bound	= YADE_PTR_CAST<Bound>(aabb);
	rootBody->physicalParameters 	= physics;
	
}

YADE_PLUGIN((ThreePointBending));

YADE_REQUIRE_FEATURE(PHYSPAR);

