/*************************************************************************
*  Copyright (C) 2007 by Bruno Chareyre		                         *
*  bruno.chareyre@hmg.inpg.fr                                            *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "CohesiveTriaxialTest.hpp"

#include<yade/pkg-dem/CohesiveFrictionalContactLaw.hpp>
#include<yade/pkg-dem/Ip2_2xCohFrictMat_CohFrictPhys.hpp>
#include<yade/pkg-dem/CohFrictMat.hpp>
#include<yade/pkg-dem/GlobalStiffnessTimeStepper.hpp>

#include<yade/pkg-dem/TriaxialStressController.hpp>
#include<yade/pkg-dem/TriaxialCompressionEngine.hpp>
#include <yade/pkg-dem/TriaxialStateRecorder.hpp>

#include<yade/pkg-common/Aabb.hpp>
#include<yade/core/Scene.hpp>
#include<yade/pkg-common/InsertionSortCollider.hpp>
#include<yade/core/Interaction.hpp>
#include<yade/pkg-common/BoundDispatcher.hpp>

#include<yade/pkg-common/GravityEngines.hpp>
#include<yade/pkg-dem/NewtonIntegrator.hpp>
#include<yade/pkg-dem/Ig2_Sphere_Sphere_ScGeom.hpp>
#include<yade/pkg-dem/Ig2_Box_Sphere_ScGeom.hpp>

#include<yade/pkg-common/InteractionGeometryDispatcher.hpp>
#include<yade/pkg-common/InteractionPhysicsDispatcher.hpp>
#include<yade/core/Body.hpp>
#include<yade/pkg-common/Box.hpp>
#include<yade/pkg-common/Sphere.hpp>

#include<yade/pkg-common/ForceResetter.hpp>

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



using namespace boost;
using namespace std;


typedef pair<Vector3r, Real> BasicSphere;
//! make a list of spheres non-overlapping sphere
string GenerateCloud_cohesive(vector<BasicSphere>& sphere_list, Vector3r lowerCorner, Vector3r upperCorner, long number, Real rad_std_dev, Real porosity);

/*
CohesiveTriaxialTest::CohesiveTriaxialTest () : FileGenerator()
{
	lowerCorner 		= Vector3r(0,0,0);
	upperCorner 		= Vector3r(1,1,1);
	thickness 		= 0.001;
	importFilename 		= "";
	outputFileName 		= "./CohesiveTriaxialTest.xml";
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
	WallStressRecordFile = "./WallStresses";

	rotationBlocked = false;
	//	boxWalls 		= false;
	boxWalls 		= true;
	internalCompaction	=false;

	dampingForce = 0.2;
	dampingMomentum = 0.2;
	defaultDt = 0.001;
	
	timeStepUpdateInterval = 50;
	timeStepOutputInterval = 50;
	wallStiffnessUpdateInterval = 10;
	radiusControlInterval = 10;
	numberOfGrains = 400;
	radiusDeviation = 0.3;
	strainRate = 10;
	StabilityCriterion = 0.01;
	autoCompressionActivation = false;
	maxMultiplier = 1.01;
	finalMaxMultiplier = 1.0001;
	
	sphereYoungModulus  = 15000000.0;
	sphereKsDivKn  = 0.5;
	sphereFrictionDeg   = 18.0;
	normalCohesion = 0;
	shearCohesion = 0;
	setCohesionOnNewContacts = false;
	density			= 2600;
	
	boxYoungModulus   = 15000000.0;
	boxKsDivKn  = 0.2;
	boxFrictionDeg   = 0.f;
	gravity 	= Vector3r(0,-9.81,0);
	
	sigma_iso = 50000;
	
//	wall_top_id =0;
// 	wall_bottom_id =0;
// 	wall_left_id =0;
// 	all_right_id =0;
// 	wall_front_id =0;
// 	wall_back_id =0;
}*/


CohesiveTriaxialTest::~CohesiveTriaxialTest ()
{

}




bool CohesiveTriaxialTest::generate()
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
//			forcerec->startId = body->getId();
//			forcerec->endId   = body->getId();
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
	else message=GenerateCloud_cohesive(sphere_list, lowerCorner, upperCorner, numberOfGrains, radiusDeviation, 0.75);
	
	vector<BasicSphere>::iterator it = sphere_list.begin();
	vector<BasicSphere>::iterator it_end = sphere_list.end();
			
	for (;it!=it_end; ++it)
	{
		cerr << "sphere (" << it->first << " " << it->second << endl;
		createSphere(body,it->first,it->second,true);
		rootBody->bodies->insert(body);
	}
	
	return true;

}


void CohesiveTriaxialTest::createSphere(shared_ptr<Body>& body, Vector3r position, Real radius, bool dynamic )
{
	body = shared_ptr<Body>(new Body(body_id_t(0),2));
	shared_ptr<CohFrictMat> physics(new CohFrictMat);
	shared_ptr<Aabb> aabb(new Aabb);
	shared_ptr<Sphere> iSphere(new Sphere);
	
	Quaternionr q(Mathr::SymmetricRandom(),Mathr::SymmetricRandom(),Mathr::SymmetricRandom(),Mathr::SymmetricRandom());
	q.normalize();
	
	body->isDynamic			= dynamic;
	
	body->state->angVel		= Vector3r(0,0,0);
	body->state->vel		= Vector3r(0,0,0);
	body->state->mass		= 4.0/3.0*Mathr::PI*radius*radius*radius*density;
	
	body->state->inertia		= Vector3r( 	2.0/5.0*body->state->mass*radius*radius,
							2.0/5.0*body->state->mass*radius*radius,
   							2.0/5.0*body->state->mass*radius*radius);
	body->state->se3			= Se3r(position,q);
	physics->young			= sphereYoungModulus;
	physics->poisson		= sphereKsDivKn;
	physics->frictionAngle		= sphereFrictionDeg * Mathr::PI/180.0;

	if((!dynamic) && (!boxWalls))
	{
		physics->young			= boxYoungModulus;
		physics->poisson		= boxKsDivKn;
		physics->frictionAngle		= boxFrictionDeg * Mathr::PI/180.0;
	}
	
	aabb->diffuseColor		= Vector3r(0,1,0);
	
	iSphere->radius			= radius;
	iSphere->color		= Vector3r(Mathr::UnitRandom(),Mathr::UnitRandom(),Mathr::UnitRandom());
	iSphere->wire			= false;
	
	body->shape	= iSphere;
	body->bound		= aabb;
	body->material	= physics;
}


void CohesiveTriaxialTest::createBox(shared_ptr<Body>& body, Vector3r position, Vector3r extents, bool wire)
{
	body = shared_ptr<Body>(new Body(body_id_t(0),2));
	shared_ptr<CohFrictMat> physics(new CohFrictMat);
	shared_ptr<Aabb> aabb(new Aabb);

	shared_ptr<Box> iBox(new Box);
	
	body->isDynamic			= false;
	
	body->state->angVel		= Vector3r(0,0,0);
	body->state->vel		= Vector3r(0,0,0);
	body->state->mass			= 0; 
	//physics->mass			= extents[0]*extents[1]*extents[2]*density*2; 
	body->state->inertia		= Vector3r(
						body->state->mass*(extents[1]*extents[1]+extents[2]*extents[2])/3
						, body->state->mass*(extents[0]*extents[0]+extents[2]*extents[2])/3
						, body->state->mass*(extents[1]*extents[1]+extents[0]*extents[0])/3
						);
//	physics->mass			= 0;
//	physics->inertia		= Vector3r(0,0,0);
	body->state->se3			= Se3r(position,Quaternionr::Identity());

	physics->young			= boxYoungModulus;
	physics->poisson		= boxKsDivKn;
	physics->frictionAngle		= boxFrictionDeg * Mathr::PI/180.0;

	aabb->diffuseColor		= Vector3r(1,0,0);
	
	iBox->extents			= extents;
	iBox->color		= Vector3r(1,1,1);
	iBox->wire			= wire;

	body->bound		= aabb;
	body->shape	= iBox;
	body->material	= physics;	
}


void CohesiveTriaxialTest::createActors(shared_ptr<Scene>& rootBody)
{
	
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

	
	shared_ptr<NewtonIntegrator> newton(new NewtonIntegrator);
	newton->damping=dampingMomentum;

	
	shared_ptr<GlobalStiffnessTimeStepper> globalStiffnessTimeStepper(new GlobalStiffnessTimeStepper);
	globalStiffnessTimeStepper->timeStepUpdateInterval = timeStepUpdateInterval;
	globalStiffnessTimeStepper->defaultDt = defaultDt;
	globalStiffnessTimeStepper->timestepSafetyCoefficient = 0.2;
	shared_ptr<CohesiveFrictionalContactLaw> cohesiveFrictionalContactLaw(new CohesiveFrictionalContactLaw);
	triaxialcompressionEngine = shared_ptr<TriaxialCompressionEngine> (new TriaxialCompressionEngine);
	triaxialcompressionEngine-> stiffnessUpdateInterval = wallStiffnessUpdateInterval;// = stiffness update interval
	triaxialcompressionEngine-> radiusControlInterval = radiusControlInterval;// = stiffness update interval
	triaxialcompressionEngine-> sigma_iso = sigmaIsoCompaction;
	triaxialcompressionEngine-> sigmaLateralConfinement = sigmaLateralConfinement;
	triaxialcompressionEngine-> sigmaIsoCompaction = sigmaIsoCompaction;
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
	triaxialStateRecorder-> file 	= WallStressRecordFile;
	triaxialStateRecorder-> iterPeriod 		= recordIntervalIter;
	//triaxialStateRecorder-> thickness 		= thickness;
	
	rootBody->engines.clear();
	rootBody->engines.push_back(shared_ptr<Engine>(new ForceResetter));
	rootBody->engines.push_back(boundDispatcher);
	rootBody->engines.push_back(shared_ptr<Engine>(new InsertionSortCollider));
	rootBody->engines.push_back(interactionGeometryDispatcher);
	rootBody->engines.push_back(interactionPhysicsDispatcher);
	rootBody->engines.push_back(cohesiveFrictionalContactLaw);
	rootBody->engines.push_back(triaxialcompressionEngine);
	rootBody->engines.push_back(globalStiffnessTimeStepper);
	rootBody->engines.push_back(triaxialStateRecorder);
	//rootBody->engines.push_back(hydraulicForceEngine);//<-------------HYDRAULIC ENGINE HERE
	rootBody->engines.push_back(newton);

	rootBody->initializers.clear();
	rootBody->initializers.push_back(boundDispatcher);
	
}


void CohesiveTriaxialTest::positionRootBody(shared_ptr<Scene>& rootBody)
{	
	shared_ptr<Aabb> aabb(new Aabb);
	aabb->diffuseColor		= Vector3r(0,0,1);
}



string GenerateCloud_cohesive(vector<BasicSphere>& sphere_list, Vector3r lowerCorner, Vector3r upperCorner, long number, Real rad_std_dev, Real porosity)
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
		
	Real mean_radius = std::pow(dimensions.x()*dimensions.y()*dimensions.z()*(1-porosity)/(3.1416*1.3333*number),0.333333);
        //cerr << mean_radius;

	std::cerr << "generating aggregates ... ";
	
	long t, i;
	for (i=0; i<number; ++i) {
		BasicSphere s;
		s.second = (random1()-0.5)*rad_std_dev*mean_radius+mean_radius;
		for (t=0; t<tries; ++t) {
			s.first.x() = lowerCorner.x()+s.second+(dimensions.x()-2*s.second)*random1();
			s.first.y() = lowerCorner.y()+s.second+(dimensions.y()-2*s.second)*random1();
			s.first.z() = lowerCorner.z()+s.second+(dimensions.z()-2*s.second)*random1();
			bool overlap=false;
			for (long j=0; (j<i && !overlap); j++)
				if ( pow(sphere_list[j].second+s.second, 2) > (sphere_list[j].first-s.first).squaredNorm()) overlap=true;
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



YADE_PLUGIN((CohesiveTriaxialTest));

