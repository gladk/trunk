/*************************************************************************
*  Copyright (C) 2006 by Bruno Chareyre		                         *
*  bruno.chareyre@hmg.inpg.fr                                            *
*  Copyright (C) 2006 by luc Scholtes                                    *
*  luc.scholtes@hmg.inpg.fr                                              *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
#include "CapillaryTriaxialTest.hpp"

#include<yade/pkg-dem/ElasticContactLaw.hpp>
#include <yade/pkg-dem/Law2_ScGeom_CapillaryPhys_Capillarity.hpp>
// #include<yade/pkg-dem/Ip2_FrictMat_FrictMat_FrictPhys.hpp>
#include<yade/pkg-dem/Ip2_FrictMat_FrictMat_CapillaryPhys.hpp>
#include<yade/pkg-common/ElastMat.hpp>
#include<yade/pkg-dem/GlobalStiffnessTimeStepper.hpp>

#include<yade/pkg-dem/TriaxialStressController.hpp>
#include<yade/pkg-dem/TriaxialCompressionEngine.hpp>
#include <yade/pkg-dem/TriaxialStateRecorder.hpp>
#include <yade/pkg-dem/CapillaryStressRecorder.hpp>

#include<yade/pkg-common/Aabb.hpp>
#include<yade/core/Scene.hpp>
#include<yade/pkg-common/InsertionSortCollider.hpp>
#include<yade/core/Interaction.hpp>
#include<yade/pkg-common/BoundDispatcher.hpp>
#include<yade/pkg-common/Bo1_Sphere_Aabb.hpp>

#include<yade/pkg-common/GravityEngines.hpp>
#include<yade/pkg-dem/NewtonIntegrator.hpp>

#include<yade/pkg-common/InteractionGeometryDispatcher.hpp>
#include<yade/pkg-common/InteractionPhysicsDispatcher.hpp>
#include<yade/core/Body.hpp>
#include<yade/pkg-common/Box.hpp>
#include<yade/pkg-common/Sphere.hpp>
#include<yade/pkg-dem/Ig2_Sphere_Sphere_ScGeom.hpp>

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
//! generate a list of non-overlapping spheres
string GenerateCloud_water(vector<BasicSphere>& sphere_list, Vector3r lowerCorner, Vector3r upperCorner, long number, Real rad_std_dev, Real porosity);


// CapillaryTriaxialTest::CapillaryTriaxialTest () : FileGenerator()
// {
// 	lowerCorner 		= Vector3r(0,0,0);
// 	upperCorner 		= Vector3r(0.001,0.001,0.001);
// 	thickness 		= 0.00001;
// 	importFilename 		= ""; // "./small.sdec.xyz";
// 	Key			="";
// 	outputFileName 		= "./CapillaryTriaxialTest"+Key+".xml";
// 	//nlayers = 1;
// 	wall_top 		= true;
// 	wall_bottom 		= true;
// 	wall_1			= true;
// 	wall_2			= true;
// 	wall_3			= true;
// 	wall_4			= true;
// 	wall_top_wire 		= true;
// 	wall_bottom_wire	= true;
// 	wall_1_wire		= true;
// 	wall_2_wire		= true;
// 	wall_3_wire		= true;
// 	wall_4_wire		= true;
// 	spheresColor		= Vector3r(0.8,0.3,0.3);
// 	spheresRandomColor	= false;
// 	recordIntervalIter	= 20;
// 	WallStressRecordFile = "./WallStresses"+Key;
// 	capillaryStressRecordFile	= "./CapillaryStresses"+Key;
// 
// 	rotationBlocked = false;
// 	//	boxWalls 		= false;
// 	boxWalls 		= true;
// 	internalCompaction	=false;
// 	
// 	dampingForce = 0.2;
// 	dampingMomentum = 0.2;
// 	defaultDt = 0.0001;
// 	
// 	timeStepUpdateInterval = 50;
// 	timeStepOutputInterval = 50;
// 	wallStiffnessUpdateInterval = 10;
// 	radiusControlInterval = 10;
// 	numberOfGrains = 400;
// 	Rdispersion = 0.667;
// 	strainRate = 10;
// 	maxWallVelocity=0.1;
// 	StabilityCriterion = 0.01;
// 	autoCompressionActivation = true;
// 	maxMultiplier = 1.01;
// 	finalMaxMultiplier = 1.001;
// 	
// 	sphereYoungModulus  = 5000000.0;
// 	sphereKsDivKn  = 0.5;	
// 	sphereFrictionDeg = 30.0;
// 	compactionFrictionDeg   = sphereFrictionDeg;
// 	density			= 2600;
// 	
// 	boxYoungModulus   = 5000000.0;
// 	boxKsDivKn  = 0.2;
// 	boxFrictionDeg   = 0.f;
// 	gravity 	= Vector3r(0,-9.81,0);
// 	
// 	sigmaIsoCompaction = 5000;
// 	sigmaLateralConfinement=sigmaIsoCompaction;
// 
// 	water = true;
// 	CapillaryPressure = 0;
//  	fusionDetection = 1;
// 	binaryFusion = 1;
// 
// 	wallOversizeFactor=1.3;
// 	
// //	wall_top_id =0;
// // 	wall_bottom_id =0;
// // 	wall_left_id =0;
// // 	all_right_id =0;
// // 	wall_front_id =0;
// // 	wall_back_id =0;
// }


CapillaryTriaxialTest::~CapillaryTriaxialTest ()
{

}



bool CapillaryTriaxialTest::generate()
{
//	unsigned int startId=boost::numeric::bounds<unsigned int>::highest(), endId=0; // record forces from group 2
	message="";	
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
	 						wallOversizeFactor*fabs(lowerCorner[0]-upperCorner[0])/2+thickness,
							thickness/2.0,
	 						wallOversizeFactor*fabs(lowerCorner[2]-upperCorner[2])/2+thickness);
	
		createBox(body,center,halfSize,wall_bottom_wire);
	 	if(wall_bottom) {
			rootBody->bodies->insert(body);
			triaxialcompressionEngine->wall_bottom_id = body->getId();
			//triaxialStateRecorder->wall_bottom_id = body->getId();
			}
		//forcerec->id = body->getId();
	
	// top box
	 	center			= Vector3r(
	 						(lowerCorner[0]+upperCorner[0])/2,
	 						upperCorner[1]+thickness/2.0,
	 						(lowerCorner[2]+upperCorner[2])/2);
	 	halfSize		= Vector3r(
	 						wallOversizeFactor*fabs(lowerCorner[0]-upperCorner[0])/2+thickness,
	 						thickness/2.0,
	 						wallOversizeFactor*fabs(lowerCorner[2]-upperCorner[2])/2+thickness);
	
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
	 						wallOversizeFactor*fabs(lowerCorner[1]-upperCorner[1])/2+thickness,
	 						wallOversizeFactor*fabs(lowerCorner[2]-upperCorner[2])/2+thickness);
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
	 						wallOversizeFactor*fabs(lowerCorner[1]-upperCorner[1])/2+thickness,
	 						wallOversizeFactor*fabs(lowerCorner[2]-upperCorner[2])/2+thickness);
	 	
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
	 						wallOversizeFactor*fabs(lowerCorner[0]-upperCorner[0])/2+thickness,
	 						wallOversizeFactor*fabs(lowerCorner[1]-upperCorner[1])/2+thickness,
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
	 						wallOversizeFactor*fabs(lowerCorner[0]-upperCorner[0])/2+thickness,
	 						wallOversizeFactor*fabs(lowerCorner[1]-upperCorner[1])/2+thickness,
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
	else message+=GenerateCloud_water(sphere_list, lowerCorner, upperCorner, numberOfGrains, Rdispersion, 0.75);
	
	vector<BasicSphere>::iterator it = sphere_list.begin();
	vector<BasicSphere>::iterator it_end = sphere_list.end();
			
	for (;it!=it_end; ++it)
	{
		cerr << "sphere (" << it->first << " " << it->second << ")"<<endl;
		createSphere(body,it->first,it->second,false,true);
		rootBody->bodies->insert(body);
	}	
	
	return true;
//  	return "Generated a sample inside box of dimensions: (" 
//  		+ lexical_cast<string>(lowerCorner[0]) + "," 
//  		+ lexical_cast<string>(lowerCorner[1]) + "," 
//  		+ lexical_cast<string>(lowerCorner[2]) + ") and (" 
//  		+ lexical_cast<string>(upperCorner[0]) + "," 
//  		+ lexical_cast<string>(upperCorner[1]) + "," 
//  		+ lexical_cast<string>(upperCorner[2]) + ").";

}


void CapillaryTriaxialTest::createSphere(shared_ptr<Body>& body, Vector3r position, Real radius, bool big, bool dynamic )
{
	body = shared_ptr<Body>(new Body(body_id_t(0),2));
	shared_ptr<FrictMat> physics(new FrictMat);
	shared_ptr<Aabb> aabb(new Aabb);

	shared_ptr<Sphere> iSphere(new Sphere);
	
	Quaternionr q(Mathr::SymmetricRandom(),Mathr::SymmetricRandom(),Mathr::SymmetricRandom(),Mathr::SymmetricRandom());
	q.normalize();
//	q.FromAxisAngle( Vector3r(0,0,1),0);
	
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


void CapillaryTriaxialTest::createBox(shared_ptr<Body>& body, Vector3r position, Vector3r extents, bool wire)
{
	body = shared_ptr<Body>(new Body(body_id_t(0),2));
	shared_ptr<FrictMat> physics(new FrictMat);
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


void CapillaryTriaxialTest::createActors(shared_ptr<Scene>& rootBody)
{
	
	Real distanceFactor = 1.3;//Create potential interactions as soon as the distance is less than factor*(rad1+rad2) 
	
	shared_ptr<InteractionGeometryDispatcher> interactionGeometryDispatcher(new InteractionGeometryDispatcher);
	
	shared_ptr<Ig2_Sphere_Sphere_ScGeom> iS2IS4SContactGeometry(new Ig2_Sphere_Sphere_ScGeom);
	iS2IS4SContactGeometry->interactionDetectionFactor = distanceFactor;//Detect potential distant interaction (meniscii)
	
	interactionGeometryDispatcher->add(iS2IS4SContactGeometry);
	interactionGeometryDispatcher->add("Ig2_Box_Sphere_ScGeom");

	shared_ptr<InteractionPhysicsDispatcher> interactionPhysicsDispatcher(new InteractionPhysicsDispatcher);
//	interactionPhysicsDispatcher->add("Ip2_FrictMat_FrictMat_FrictPhys");
// Unhandled exception: St13runtime_error : Class `Ip2_FrictMat_FrictMat_FrictPhys' could not be cast to required 2D Functor	
	/// OLD
	//interactionPhysicsDispatcher->add("BodyMacroParameters","BodyMacroParameters","MacroMicroElasticRelationshipsWater");
	/// NEW
	shared_ptr<InteractionPhysicsFunctor> ss(new Ip2_FrictMat_FrictMat_CapillaryPhys);
	interactionPhysicsDispatcher->add(ss);
	
		
	shared_ptr<BoundDispatcher> boundDispatcher	= shared_ptr<BoundDispatcher>(new BoundDispatcher);
	
	shared_ptr<Bo1_Sphere_Aabb> interactingSphere2AABB(new Bo1_Sphere_Aabb);
	interactingSphere2AABB->aabbEnlargeFactor = distanceFactor;//Detect potential distant interaction (meniscii)
	
	boundDispatcher->add(interactingSphere2AABB);
	boundDispatcher->add("Bo1_Box_Aabb");

	shared_ptr<GravityEngine> gravityCondition(new GravityEngine);
	gravityCondition->gravity = gravity;
	
	shared_ptr<GlobalStiffnessTimeStepper> globalStiffnessTimeStepper(new GlobalStiffnessTimeStepper);
	globalStiffnessTimeStepper->timeStepUpdateInterval = timeStepUpdateInterval;
	globalStiffnessTimeStepper->defaultDt = defaultDt;
	
	shared_ptr<ElasticContactLaw> elasticContactLaw(new ElasticContactLaw);
	elasticContactLaw->neverErase = true;
	

	// capillary
	shared_ptr<Law2_ScGeom_CapillaryPhys_Capillarity> capillaryCohesiveLaw(new Law2_ScGeom_CapillaryPhys_Capillarity); 
	capillaryCohesiveLaw->CapillaryPressure = CapillaryPressure;
// 	capillaryCohesiveLaw->fusionDetection = fusionDetection;
// 	capillaryCohesiveLaw->binaryFusion = binaryFusion;

	//shared_ptr<StiffnessCounter> stiffnesscounter(new StiffnessCounter);
	//stiffnesscounter->interval = timeStepUpdateInterval;
	
	
	// moving walls to regulate the stress applied + compress when the packing is dense an stable
	//cerr << "triaxialcompressionEngine = shared_ptr<TriaxialCompressionEngine> (new TriaxialCompressionEngine);" << std::endl;
	triaxialcompressionEngine = shared_ptr<TriaxialCompressionEngine> (new TriaxialCompressionEngine);
	triaxialcompressionEngine-> stiffnessUpdateInterval = wallStiffnessUpdateInterval;// = stiffness update interval
	triaxialcompressionEngine-> radiusControlInterval = radiusControlInterval;// = stiffness update interval
	triaxialcompressionEngine-> sigmaIsoCompaction = sigmaIsoCompaction;
	triaxialcompressionEngine-> sigmaLateralConfinement = sigmaLateralConfinement;
	triaxialcompressionEngine->max_vel = maxWallVelocity;
	triaxialcompressionEngine-> thickness = thickness;
	triaxialcompressionEngine->strainRate = strainRate;
	triaxialcompressionEngine->StabilityCriterion = StabilityCriterion;
	triaxialcompressionEngine->autoCompressionActivation = autoCompressionActivation;
	triaxialcompressionEngine->internalCompaction = internalCompaction;
	triaxialcompressionEngine->maxMultiplier = maxMultiplier;
	triaxialcompressionEngine->finalMaxMultiplier = finalMaxMultiplier;
	triaxialcompressionEngine->Key = Key;
	triaxialcompressionEngine->frictionAngleDegree = sphereFrictionDeg;
		
	//cerr << "fin de section triaxialcompressionEngine = shared_ptr<TriaxialCompressionEngine> (new TriaxialCompressionEngine);" << std::endl;
	
	// recording global stress
	triaxialStateRecorder = shared_ptr<TriaxialStateRecorder>(new TriaxialStateRecorder);
	triaxialStateRecorder-> file 		= WallStressRecordFile + Key;
	triaxialStateRecorder-> iterPeriod	= recordIntervalIter;
	//triaxialStateRecorderer-> thickness 		= thickness;

	// recording capillary stress
	capillaryStressRecorder = shared_ptr<CapillaryStressRecorder>(new CapillaryStressRecorder);
	capillaryStressRecorder -> file		= capillaryStressRecordFile + Key;
	capillaryStressRecorder -> iterPeriod	= recordIntervalIter;

	#if 0	
	// moving walls to regulate the stress applied
	//cerr << "triaxialstressController = shared_ptr<TriaxialStressController> (new TriaxialStressController);" << std::endl;
	triaxialstressController = shared_ptr<TriaxialStressController> (new TriaxialStressController);
	triaxialstressController-> stiffnessUpdateInterval = 20;// = recordIntervalIter
	triaxialstressController-> sigma_iso = sigma_iso;
	triaxialstressController-> max_vel = 0.0001;
	triaxialstressController-> thickness = thickness;
	triaxialstressController->wall_bottom_activated = false;
	triaxialstressController->wall_top_activated = false;	
		//cerr << "fin de sezction triaxialstressController = shared_ptr<TriaxialStressController> (new TriaxialStressController);" << std::endl;
	#endif
	
	
	rootBody->engines.clear();
	rootBody->engines.push_back(shared_ptr<Engine>(new ForceResetter));
//	rootBody->engines.push_back(sdecTimeStepper);	
	rootBody->engines.push_back(boundDispatcher);
	rootBody->engines.push_back(shared_ptr<Engine>(new InsertionSortCollider));
	rootBody->engines.push_back(interactionGeometryDispatcher);
	rootBody->engines.push_back(interactionPhysicsDispatcher);
	rootBody->engines.push_back(elasticContactLaw);

	// capillary
	
	
	//rootBody->engines.push_back(stiffnesscounter);
	//rootBody->engines.push_back(stiffnessMatrixTimeStepper);
	rootBody->engines.push_back(globalStiffnessTimeStepper);
	if(water)
	{
		rootBody->engines.push_back(capillaryCohesiveLaw);
		rootBody->engines.push_back(capillaryStressRecorder);
	}
	rootBody->engines.push_back(triaxialcompressionEngine);
	rootBody->engines.push_back(triaxialStateRecorder);
	//rootBody->engines.push_back(gravityCondition);
	
	rootBody->engines.push_back(shared_ptr<Engine> (new NewtonIntegrator));
	
	//if(!rotationBlocked)
	//	rootBody->engines.push_back(orientationIntegrator);
	//rootBody->engines.push_back(triaxialstressController);
	
	rootBody->initializers.clear();
	rootBody->initializers.push_back(boundDispatcher);
	
}


void CapillaryTriaxialTest::positionRootBody(shared_ptr<Scene>& rootBody)
{	
}


string GenerateCloud_water(vector<BasicSphere>& sphere_list, Vector3r lowerCorner, Vector3r upperCorner, long number, Real rad_std_dev, Real porosity)
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
        Real Rmin=mean_radius, Rmax=mean_radius;

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
			if (!overlap){
				sphere_list.push_back(s);
				Rmin = std::min(Rmin,s.second);
				Rmax = std::max(Rmax,s.second);
				break;}
		}
		if (t==tries) return "More than " + lexical_cast<string>(tries) +
					" tries while generating sphere number " +
					lexical_cast<string>(i+1) + "/" + lexical_cast<string>(number) + ".";
	}
	return "Generated a sample with " + lexical_cast<string>(number) + "spheres inside box of dimensions: (" 
			+ lexical_cast<string>(dimensions[0]) + "," 
			+ lexical_cast<string>(dimensions[1]) + "," 
			+ lexical_cast<string>(dimensions[2]) + ")."
			+ "  mean radius=" + lexical_cast<string>(mean_radius) +
			+ "  Rmin =" + lexical_cast<string>(Rmin) +
			+ "  Rmax =" + lexical_cast<string>(Rmax) + ".";
}



YADE_PLUGIN((CapillaryTriaxialTest));

//YADE_REQUIRE_FEATURE(PHYSPAR);

