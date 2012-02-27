/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*  Copyright (C) 2008 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*  Copyright (C) 2007 by Bruno Chareyre		                         *
*  bruno.chareyre@hmg.inpg.fr                                            *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "SnowCreepTest.hpp"

#include<yade/pkg-dem/CohesiveFrictionalContactLaw.hpp>
#include<yade/pkg-dem/Ip2_2xCohFrictMat_CohFrictPhys.hpp>
#include<yade/pkg-dem/CohFrictMat.hpp>
#include<yade/pkg-dem/SDECLinkPhysics.hpp>
#include<yade/pkg-dem/GlobalStiffnessTimeStepper.hpp>
#include<yade/pkg-dem/PositionOrientationRecorder.hpp>

#include<yade/pkg-dem/AveragePositionRecorder.hpp>
#include<yade/pkg-dem/ForceRecorder.hpp>
#include<yade/pkg-dem/VelocityRecorder.hpp>
#include<yade/pkg-dem/TriaxialStressController.hpp>
#include<yade/pkg-dem/TriaxialCompressionEngine.hpp>
#include <yade/pkg-dem/TriaxialStateRecorder.hpp>

#include<yade/pkg-common/BoxModel.hpp>
#include<yade/pkg-common/Aabb.hpp>
#include<yade/pkg-common/SphereModel.hpp>
#include<yade/core/Scene.hpp>
#include<yade/pkg-common/InsertionSortCollider.hpp>
#include<yade/lib-serialization/IOFormatManager.hpp>
#include<yade/core/Interaction.hpp>
#include<yade/pkg-common/BoundDispatcher.hpp>

#include<yade/pkg-common/GravityEngines.hpp>
#include<yade/pkg-dem/Ig2_Sphere_Sphere_ScGeom.hpp>
#include<yade/pkg-dem/Ig2_Box_Sphere_ScGeom.hpp>
#include<yade/pkg-common/PhysicalActionApplier.hpp>
#include<yade/pkg-common/PhysicalActionDamper.hpp>
#include<yade/pkg-common/CundallNonViscousDamping.hpp>
#include<yade/pkg-common/CundallNonViscousDamping.hpp>

#include<yade/pkg-common/InteractionGeometryDispatcher.hpp>
#include<yade/pkg-common/InteractionPhysicsDispatcher.hpp>
#include<yade/core/Body.hpp>
#include<yade/pkg-common/Box.hpp>
#include<yade/pkg-common/Sphere.hpp>

#include<yade/pkg-common/ForceResetter.hpp>

#include<yade/pkg-common/StateMetaEngine.hpp>


#include<yade/pkg-dem/Shop.hpp>

#include <boost/filesystem/convenience.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/numeric/conversion/bounds.hpp>
#include <boost/limits.hpp>

// random
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>

YADE_REQUIRE_FEATURE(geometricalmodel)

using namespace boost;
using namespace std;


typedef pair<Vector3r, Real> BasicSphere;
//! make a list of spheres non-overlapping sphere
string GenerateCloud_snow(vector<BasicSphere>& sphere_list, Vector3r lowerCorner, Vector3r upperCorner, long number, Real rad_std_dev, Real porosity);


SnowCreepTest::SnowCreepTest () : FileGenerator()
{
	lowerCorner 		= Vector3r(0,0,0);
	upperCorner 		= Vector3r(0.01,0.01,0.01);
	thickness 		= 0.001;
	importFilename 		= "";
	outputFileName 		= "./SnowCreepTest.xml";
	//nlayers = 1;
	wall_top 		= true;
	wall_bottom 		= true;
	wall_1			= true;
	wall_2			= true;
	wall_3			= true;
	wall_4			= true;
	wall_top_wire 		= true;
	wall_bottom_wire	= true;
	wall_1_wire		= true;
	wall_2_wire		= true;
	wall_3_wire		= true;
	wall_4_wire		= true;
	spheresColor		= Vector3r(0.8,0.3,0.3);
	spheresRandomColor	= false;
	recordBottomForce	= true;
	forceRecordFile		= "./force";
	recordAveragePositions	= true;
	positionRecordFile	= "./positions";
	recordIntervalIter	= 20;
	velocityRecordFile 	= "./velocities";
	saveAnimationSnapshots = false;
	AnimationSnapshotsBaseName = "./snapshots/snap";
	WallStressRecordFile = "./WallStresses";

	rotationBlocked = false;
	//	boxWalls 		= false;
	boxWalls 		= true;
	internalCompaction	=false;

	dampingForce = 0.2;
	dampingMomentum = 0.2;
	defaultDt = 0.00001;
	
	timeStepUpdateInterval = 50;
	timeStepOutputInterval = 50;
	wallStiffnessUpdateInterval = 10;
	radiusControlInterval = 10;
	numberOfGrains = 100;
	radiusDeviation = 0.3;
	strainRate = 10;
	StabilityCriterion = 0.01;
	autoCompressionActivation = false;
	maxMultiplier = 1.01;
	finalMaxMultiplier = 1.0001;
	
	sphereYoungModulus  = 15000000.0;
	spherePoissonRatio  = 0.5;
	sphereFrictionDeg   = 18.0;
	normalCohesion = 500000;
	shearCohesion = 500000;
	setCohesionOnNewContacts = false;
	density			= 2600;
	
	boxYoungModulus   = 15000000.0;
	boxPoissonRatio  = 0.2;
	boxFrictionDeg   = 0.f;
	gravity 	= Vector3r(0,-9.81,0);
	
	sigma_iso = 50000;
	
//	wall_top_id =0;
// 	wall_bottom_id =0;
// 	wall_left_id =0;
// 	all_right_id =0;
// 	wall_front_id =0;
// 	wall_back_id =0;


	creep_viscosity = 4000000;
}


SnowCreepTest::~SnowCreepTest ()
{

}




bool SnowCreepTest::generate()
{
//	unsigned int startId=boost::numeric::bounds<unsigned int>::highest(), endId=0; // record forces from group 2
	
	rootBody = shared_ptr<Scene>(new Scene);
	createActors(rootBody);
	positionRootBody(rootBody);

	shared_ptr<Body> body;
	
	if(boxWalls)
	{
	// bottom box
	 	Vector3r center		= Vector3r(
	 						(lowerCorner[0]+upperCorner[0])/2,
	 						lowerCorner[1]-thickness/2.0,
	 						(lowerCorner[2]+upperCorner[2])/2);
	 	Vector3r halfSize	= Vector3r(
	 						1.5*fabs(lowerCorner[0]-upperCorner[0])/2+thickness,
							thickness/2.0,
	 						1.5*fabs(lowerCorner[2]-upperCorner[2])/2+thickness);
	
		createBox(body,center,halfSize,wall_bottom_wire);
	 	if(wall_bottom) {
			rootBody->bodies->insert(body);
			triaxialcompressionEngine->wall_bottom_id = body->getId();
			//triaxialStateRecorder->wall_bottom_id = body->getId();
			forcerec->startId = body->getId();
			forcerec->endId   = body->getId();
			}
		//forcerec->id = body->getId();
	
	// top box
	 	center			= Vector3r(
	 						(lowerCorner[0]+upperCorner[0])/2,
	 						upperCorner[1]+thickness/2.0,
	 						(lowerCorner[2]+upperCorner[2])/2);
	 	halfSize		= Vector3r(
	 						1.5*fabs(lowerCorner[0]-upperCorner[0])/2+thickness,
	 						thickness/2.0,
	 						1.5*fabs(lowerCorner[2]-upperCorner[2])/2+thickness);
	
		createBox(body,center,halfSize,wall_top_wire);
	 	if(wall_top) {
			rootBody->bodies->insert(body);
			triaxialcompressionEngine->wall_top_id = body->getId();
			//triaxialStateRecorder->wall_top_id = body->getId();
			}
	// box 1
	
	 	center			= Vector3r(
	 						lowerCorner[0]-thickness/2.0,
	 						(lowerCorner[1]+upperCorner[1])/2,
	 						(lowerCorner[2]+upperCorner[2])/2);
		halfSize		= Vector3r(
							thickness/2.0,
	 						1.5*fabs(lowerCorner[1]-upperCorner[1])/2+thickness,
	 						1.5*fabs(lowerCorner[2]-upperCorner[2])/2+thickness);
		createBox(body,center,halfSize,wall_1_wire);
	 	if(wall_1) {
			rootBody->bodies->insert(body);
			triaxialcompressionEngine->wall_left_id = body->getId();
			//triaxialStateRecorder->wall_left_id = body->getId();
			}
	// box 2
	 	center			= Vector3r(
	 						upperCorner[0]+thickness/2.0,
	 						(lowerCorner[1]+upperCorner[1])/2,
							(lowerCorner[2]+upperCorner[2])/2);
	 	halfSize		= Vector3r(
	 						thickness/2.0,
	 						1.5*fabs(lowerCorner[1]-upperCorner[1])/2+thickness,
	 						1.5*fabs(lowerCorner[2]-upperCorner[2])/2+thickness);
	 	
		createBox(body,center,halfSize,wall_2_wire);
	 	if(wall_2) {
			rootBody->bodies->insert(body);
			triaxialcompressionEngine->wall_right_id = body->getId();
			//triaxialStateRecorder->wall_right_id = body->getId();
			}
	// box 3
	 	center			= Vector3r(
	 						(lowerCorner[0]+upperCorner[0])/2,
	 						(lowerCorner[1]+upperCorner[1])/2,
	 						lowerCorner[2]-thickness/2.0);
	 	halfSize		= Vector3r(
	 						1.5*fabs(lowerCorner[0]-upperCorner[0])/2+thickness,
	 						1.5*fabs(lowerCorner[1]-upperCorner[1])/2+thickness,
	 						thickness/2.0);
		createBox(body,center,halfSize,wall_3_wire);
	 	if(wall_3) {
			rootBody->bodies->insert(body);
			triaxialcompressionEngine->wall_back_id = body->getId();
			//triaxialStateRecorder->wall_back_id = body->getId();
			}
	
	// box 4
	 	center			= Vector3r(
	 						(lowerCorner[0]+upperCorner[0])/2,
	 						(lowerCorner[1]+upperCorner[1])/2,
	 						upperCorner[2]+thickness/2.0);
	 	halfSize		= Vector3r(
	 						1.5*fabs(lowerCorner[0]-upperCorner[0])/2+thickness,
	 						1.5*fabs(lowerCorner[1]-upperCorner[1])/2+thickness,
	 						thickness/2.0);
		createBox(body,center,halfSize,wall_3_wire);
	 	if(wall_4) {
			rootBody->bodies->insert(body);
			triaxialcompressionEngine->wall_front_id = body->getId();
			//triaxialStateRecorder->wall_front_id = body->getId();
			}
			 
	}
	
	vector<BasicSphere> sphere_list;
	if(importFilename!="") sphere_list=Shop::loadSpheresFromFile(importFilename,lowerCorner,upperCorner);
	else message=GenerateCloud_snow(sphere_list, lowerCorner, upperCorner, numberOfGrains, radiusDeviation, 0.75);
	
	vector<BasicSphere>::iterator it = sphere_list.begin();
	vector<BasicSphere>::iterator it_end = sphere_list.end();
			
	for (;it!=it_end; ++it)
	{
		cerr << "sphere (" << it->first << " " << it->second << endl;
		createSphere(body,it->first,it->second,true);
		rootBody->bodies->insert(body);
	}
	
// 	if(importFilename.size() != 0 && filesystem::exists(importFilename) )
// 	{
// 		
// 		Vector3r layersDistance (Vector3r::Zero()); 
// 		for (int layer=1; layer <= nlayers; ++layer)
// 		{			
// 			ifstream loadFile(importFilename.c_str());
// 			long int i=0;
// 			Real f,g,x,y,z,radius;
// 			while( ! loadFile.eof() )
// 			{
// 				++i;
// 				loadFile >> x;
// 				loadFile >> y;
// 				loadFile >> z;
// 				Vector3r position = (Vector3r(x,z,y) + layersDistance);
// 				loadFile >> radius;
// 			
// 				loadFile >> f;
// 				loadFile >> g;
// 				if( boxWalls ? f>1 : false ) // skip loading of SDEC walls
// 					continue;
// 				if(f==8)
// 					continue;
// 	
// 		//		if( i % 100 == 0 ) // FIXME - should display a progress BAR !!
// 		//			cout << "loaded: " << i << endl;
// 				if(f==1)
// 				{
// 					lowerCorner[0] = min(position[0]-radius , lowerCorner[0]);
// 					lowerCorner[1] = min(position[1]-radius , lowerCorner[1]);
// 					lowerCorner[2] = min(position[2]-radius , lowerCorner[2]);
// 					upperCorner[0] = max(position[0]+radius , upperCorner[0]);
// 					upperCorner[1] = max(position[1]+radius , upperCorner[1]);
// 					upperCorner[2] = max(position[2]+radius , upperCorner[2]);
// 				}
// 				createSphere(body,position,radius,false,f==1);
// 				rootBody->bodies->insert(body);
// 				if(f == 2)
// 				{
// 					startId = std::min(body->getId() , startId);
// 					endId   = std::max(body->getId() , endId);
// 				}
// 					
// 			}
// 			layersDistance.y() = upperCorner.y();
// 		}
// 	}

// create bigBall
	//Vector3r position = (upperCorner+lowerCorner)*0.5 + Vector3r(0,bigBallDropHeight,0);
	//createSphere(body,position,bigBallRadius,true,false);	
	//int bigId = 0;
// 	if(bigBall)
// 		rootBody->bodies->insert(body);
// 	bigId = body->getId();
	//forcerec->startId = startId;
	//forcerec->endId   = endId;
	//averagePositionRecorder->bigBallId = bigId;
	//velocityRecorder->bigBallId = bigId;

	
	
	return true;
//  	return "Generated a sample inside box of dimensions: (" 
//  		+ lexical_cast<string>(lowerCorner[0]) + "," 
//  		+ lexical_cast<string>(lowerCorner[1]) + "," 
//  		+ lexical_cast<string>(lowerCorner[2]) + ") and (" 
//  		+ lexical_cast<string>(upperCorner[0]) + "," 
//  		+ lexical_cast<string>(upperCorner[1]) + "," 
//  		+ lexical_cast<string>(upperCorner[2]) + ").";

}


void SnowCreepTest::createSphere(shared_ptr<Body>& body, Vector3r position, Real radius, bool dynamic )
{
	body = shared_ptr<Body>(new Body(body_id_t(0),2));
	shared_ptr<CohFrictMat> physics(new CohFrictMat);
	shared_ptr<Aabb> aabb(new Aabb);
	shared_ptr<SphereModel> gSphere(new SphereModel);
	shared_ptr<Sphere> iSphere(new Sphere);
	
	Quaternionr q(Mathr::SymmetricRandom(),Mathr::SymmetricRandom(),Mathr::SymmetricRandom(),Mathr::SymmetricRandom());
	q.Normalize();
//	q.FromAxisAngle( Vector3r(0,0,1),0);
	
	body->isDynamic			= dynamic;
	
	physics->angularVelocity	= Vector3r(0,0,0);
	physics->velocity		= Vector3r(0,0,0);
	physics->mass			= 4.0/3.0*Mathr::PI*radius*radius*radius*density;
	
	physics->inertia		= Vector3r( 	2.0/5.0*physics->mass*radius*radius,
							2.0/5.0*physics->mass*radius*radius,
							2.0/5.0*physics->mass*radius*radius);
	physics->se3			= Se3r(position,q);
	physics->young			= sphereYoungModulus;
	physics->poisson		= spherePoissonRatio;
	physics->frictionAngle		= sphereFrictionDeg * Mathr::PI/180.0;

	if((!dynamic) && (!boxWalls))
	{
		physics->young			= boxYoungModulus;
		physics->poisson		= boxPoissonRatio;
		physics->frictionAngle		= boxFrictionDeg * Mathr::PI/180.0;
	}
	
	aabb->diffuseColor		= Vector3r(0,1,0);


	gSphere->radius			= radius;
//	gSphere->diffuseColor		= ((int)(position[0]*400.0))%2?Vector3r(0.7,0.7,0.7):Vector3r(0.45,0.45,0.45);
	gSphere->diffuseColor		= spheresColor;
	gSphere->wire			= false;
	gSphere->shadowCaster		= true;
	
	iSphere->radius			= radius;
	iSphere->diffuseColor		= Vector3r(Mathr::UnitRandom(),Mathr::UnitRandom(),Mathr::UnitRandom());

	body->shape	= iSphere;
	body->geometricalModel		= gSphere;
	body->bound		= aabb;
	body->physicalParameters	= physics;
}


void SnowCreepTest::createBox(shared_ptr<Body>& body, Vector3r position, Vector3r extents, bool wire)
{
	body = shared_ptr<Body>(new Body(body_id_t(0),2));
	shared_ptr<CohFrictMat> physics(new CohFrictMat);
	shared_ptr<Aabb> aabb(new Aabb);
	shared_ptr<BoxModel> gBox(new BoxModel);
	shared_ptr<Box> iBox(new Box);
	
	Quaternionr q;
	q.FromAxisAngle( Vector3r(0,0,1),0);

	body->isDynamic			= false;
	
	physics->angularVelocity	= Vector3r(0,0,0);
	physics->velocity		= Vector3r(0,0,0);
	physics->mass			= 0; 
	//physics->mass			= extents[0]*extents[1]*extents[2]*density*2; 
	physics->inertia		= Vector3r(
							  physics->mass*(extents[1]*extents[1]+extents[2]*extents[2])/3
							, physics->mass*(extents[0]*extents[0]+extents[2]*extents[2])/3
							, physics->mass*(extents[1]*extents[1]+extents[0]*extents[0])/3
						);
//	physics->mass			= 0;
//	physics->inertia		= Vector3r(0,0,0);
	physics->se3			= Se3r(position,q);

	physics->young			= boxYoungModulus;
	physics->poisson		= boxPoissonRatio;
	physics->frictionAngle		= boxFrictionDeg * Mathr::PI/180.0;
	physics->isCohesive		= false;

	aabb->diffuseColor		= Vector3r(1,0,0);

	gBox->extents			= extents;
	gBox->diffuseColor		= Vector3r(1,1,1);
	gBox->wire			= wire;
	gBox->shadowCaster		= false;
	
	iBox->extents			= extents;
	iBox->diffuseColor		= Vector3r(1,1,1);

	body->bound		= aabb;
	body->shape	= iBox;
	body->geometricalModel		= gBox;
	body->physicalParameters	= physics;
}


void SnowCreepTest::createActors(shared_ptr<Scene>& rootBody)
{
// recording average positions
	averagePositionRecorder = shared_ptr<AveragePositionRecorder>(new AveragePositionRecorder);
	averagePositionRecorder -> outputFile 		= positionRecordFile;
	averagePositionRecorder -> interval 		= recordIntervalIter;
// recording forces
	forcerec = shared_ptr<ForceRecorder>(new ForceRecorder);
	forcerec -> outputFile 	= forceRecordFile;
	forcerec -> interval 	= recordIntervalIter;
// recording velocities
	velocityRecorder = shared_ptr<VelocityRecorder>(new VelocityRecorder);
	velocityRecorder-> outputFile 	= velocityRecordFile;
	velocityRecorder-> interval 	= recordIntervalIter;

	
	shared_ptr<InteractionGeometryDispatcher> interactionGeometryDispatcher(new InteractionGeometryDispatcher);
	shared_ptr<InteractionGeometryFunctor> s1(new Ig2_Sphere_Sphere_ScGeom);
	interactionGeometryDispatcher->add(s1);
	shared_ptr<InteractionGeometryFunctor> s2(new Ig2_Box_Sphere_ScGeom);
	interactionGeometryDispatcher->add(s2);

	shared_ptr<Ip2_2xCohFrictMat_CohFrictPhys> cohesiveFrictionalRelationships = shared_ptr<Ip2_2xCohFrictMat_CohFrictPhys> (new Ip2_2xCohFrictMat_CohFrictPhys);
	cohesiveFrictionalRelationships->shearCohesion = shearCohesion;
	cohesiveFrictionalRelationships->normalCohesion = normalCohesion;
	cohesiveFrictionalRelationships->setCohesionOnNewContacts = setCohesionOnNewContacts;
	shared_ptr<InteractionPhysicsDispatcher> interactionPhysicsDispatcher(new InteractionPhysicsDispatcher);
	interactionPhysicsDispatcher->add(cohesiveFrictionalRelationships);
		
	shared_ptr<BoundDispatcher> boundDispatcher	= shared_ptr<BoundDispatcher>(new BoundDispatcher);
	boundDispatcher->add("Bo1_Sphere_Aabb");
	boundDispatcher->add("Bo1_Box_Aabb");

	

		
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

	//shared_ptr<ElasticCriterionTimeStepper> sdecTimeStepper(new ElasticCriterionTimeStepper);
	//sdecTimeStepper->sdecGroupMask = 2;
	//sdecTimeStepper->timeStepUpdateInterval = timeStepUpdateInterval;
	
	//shared_ptr<StiffnessMatrixTimeStepper> stiffnessMatrixTimeStepper(new StiffnessMatrixTimeStepper);
	//stiffnessMatrixTimeStepper->sdecGroupMask = 2;
	//stiffnessMatrixTimeStepper->timeStepUpdateInterval = timeStepUpdateInterval;
	
	shared_ptr<GlobalStiffnessTimeStepper> globalStiffnessTimeStepper(new GlobalStiffnessTimeStepper);
	globalStiffnessTimeStepper->sdecGroupMask = 2;
	globalStiffnessTimeStepper->timeStepUpdateInterval = timeStepUpdateInterval;
	globalStiffnessTimeStepper->defaultDt = defaultDt;
	globalStiffnessTimeStepper->timestepSafetyCoefficient = 0.2;
	
	shared_ptr<CohesiveFrictionalContactLaw> cohesiveFrictionalContactLaw(new CohesiveFrictionalContactLaw);
	cohesiveFrictionalContactLaw->sdecGroupMask = 2;
	cohesiveFrictionalContactLaw->shear_creep = false;
	cohesiveFrictionalContactLaw->twist_creep = false;
	cohesiveFrictionalContactLaw->creep_viscosity = creep_viscosity;
	
	//shared_ptr<StiffnessCounter> stiffnesscounter(new StiffnessCounter);
	//stiffnesscounter->sdecGroupMask = 2;
	//stiffnesscounter->interval = timeStepUpdateInterval;
	
	
	// moving walls to regulate the stress applied + compress when the packing is dense an stable
	//cerr << "triaxialcompressionEngine = shared_ptr<TriaxialCompressionEngine> (new TriaxialCompressionEngine);" << std::endl;
	triaxialcompressionEngine = shared_ptr<TriaxialCompressionEngine> (new TriaxialCompressionEngine);
	triaxialcompressionEngine-> stiffnessUpdateInterval = wallStiffnessUpdateInterval;// = stiffness update interval
	triaxialcompressionEngine-> radiusControlInterval = radiusControlInterval;// = stiffness update interval
	triaxialcompressionEngine-> sigma_iso = sigma_iso;
	triaxialcompressionEngine-> sigmaLateralConfinement = sigma_iso;
	triaxialcompressionEngine-> sigmaIsoCompaction = sigma_iso;
	triaxialcompressionEngine-> max_vel = 1;
	triaxialcompressionEngine-> thickness = thickness;
	triaxialcompressionEngine->strainRate = strainRate;
	triaxialcompressionEngine->StabilityCriterion = StabilityCriterion;
	triaxialcompressionEngine->autoCompressionActivation = autoCompressionActivation;
	triaxialcompressionEngine->internalCompaction = internalCompaction;
	triaxialcompressionEngine->maxMultiplier = maxMultiplier;
	

	// recording global stress
	triaxialStateRecorder = shared_ptr<TriaxialStateRecorder>(new
	TriaxialStateRecorder);
	triaxialStateRecorder-> outputFile 	= WallStressRecordFile;
	triaxialStateRecorder-> interval 		= recordIntervalIter;
	//triaxialStateRecorder-> thickness 		= thickness;
	
	
	// moving walls to regulate the stress applied
	triaxialstressController = shared_ptr<TriaxialStressController> (new TriaxialStressController);
	triaxialstressController-> stiffnessUpdateInterval = 20;// = recordIntervalIter
	triaxialstressController-> sigma_iso = sigma_iso;
	triaxialstressController-> max_vel = 0.0001;
	triaxialstressController-> thickness = thickness;
	triaxialstressController->wall_bottom_activated = false;
	triaxialstressController->wall_top_activated = false;	
		//cerr << "fin de sezction triaxialstressController = shared_ptr<TriaxialStressController> (new TriaxialStressController);" << std::endl;
	
	rootBody->engines.clear();
	rootBody->engines.push_back(shared_ptr<Engine>(new ForceResetter));
//	rootBody->engines.push_back(sdecTimeStepper);	
	rootBody->engines.push_back(boundDispatcher);
	rootBody->engines.push_back(shared_ptr<Engine>(new InsertionSortCollider));
	rootBody->engines.push_back(interactionGeometryDispatcher);
	rootBody->engines.push_back(interactionPhysicsDispatcher);
	rootBody->engines.push_back(cohesiveFrictionalContactLaw);
	rootBody->engines.push_back(triaxialcompressionEngine);
	//rootBody->engines.push_back(stiffnesscounter);
	//rootBody->engines.push_back(stiffnessMatrixTimeStepper);
	rootBody->engines.push_back(globalStiffnessTimeStepper);
	rootBody->engines.push_back(triaxialStateRecorder);
//	rootBody->engines.push_back(hydraulicForceEngine);//<-------------HYDRAULIC ENGINE HERE
	rootBody->engines.push_back(actionDampingDispatcher);
	rootBody->engines.push_back(applyActionDispatcher);
	rootBody->engines.push_back(positionIntegrator);
	if(!rotationBlocked)
		rootBody->engines.push_back(orientationIntegrator);
	//rootBody->engines.push_back(triaxialstressController);
	
		
	//rootBody->engines.push_back(averagePositionRecorder);
	//rootBody->engines.push_back(velocityRecorder);
	//rootBody->engines.push_back(forcerec);
	
	if (saveAnimationSnapshots) {
	shared_ptr<PositionOrientationRecorder> positionOrientationRecorder(new PositionOrientationRecorder);
	positionOrientationRecorder->outputFile = AnimationSnapshotsBaseName;
	rootBody->engines.push_back(positionOrientationRecorder);}
	
	rootBody->initializers.clear();
	rootBody->initializers.push_back(boundDispatcher);
	
}


void SnowCreepTest::positionRootBody(shared_ptr<Scene>& rootBody)
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



string GenerateCloud_snow(vector<BasicSphere>& sphere_list, Vector3r lowerCorner, Vector3r upperCorner, long number, Real rad_std_dev, Real porosity)
{
	typedef boost::minstd_rand StdGenerator;
	static StdGenerator generator;
	static boost::variate_generator<StdGenerator&, boost::uniform_real<> >
			random1(generator, boost::uniform_real<>(0,1));
        //         static boost::variate_generator<StdGenerator&, boost::normal_distribution<> >
        //         randomN(generator, boost::normal_distribution<>(aggregateMeanRadius,aggregateSigmaRadius));

	sphere_list.clear();
	long tries = 1000; //nb of tries for positionning the next sphere
	Vector3r dimensions = upperCorner - lowerCorner;
		
	Real mean_radius = std::pow(dimensions.X()*dimensions.Y()*dimensions.Z()*(1-porosity)/(3.1416*1.3333*number),0.333333);
        //cerr << mean_radius;

	std::cerr << "generating aggregates ... ";
	
	long t, i;
	for (i=0; i<number; ++i) {
		BasicSphere s;
		for (t=0; t<tries; ++t) {
			s.second = (random1()-0.5)*rad_std_dev*mean_radius+mean_radius;
			s.first.X() = lowerCorner.X()+s.second+(dimensions.X()-2*s.second)*random1();
			s.first.Y() = lowerCorner.Y()+s.second+(dimensions.Y()-2*s.second)*random1();
			s.first.Z() = lowerCorner.Z()+s.second+(dimensions.Z()-2*s.second)*random1();
			bool overlap=false;
			for (long j=0; (j<i && !overlap); j++)
				if ( pow(sphere_list[j].second+s.second, 2) > (sphere_list[j].first-s.first).SquaredLength()) overlap=true;
			if (!overlap)
			{
				sphere_list.push_back(s);
				break;
			}			
		}
		if (t==tries) return "More than " + lexical_cast<string>(tries) +
					" tries while generating sphere number " +
					lexical_cast<string>(i+1) + "/" + lexical_cast<string>(number) + ".";
	}
	return "Generated a sample with " + lexical_cast<string>(number) + "spheres inside box of dimensions: (" 
			+ lexical_cast<string>(dimensions[0]) + "," 
			+ lexical_cast<string>(dimensions[1]) + "," 
			+ lexical_cast<string>(dimensions[2]) + ").";
}



YADE_PLUGIN((SnowCreepTest));

YADE_REQUIRE_FEATURE(PHYSPAR);

