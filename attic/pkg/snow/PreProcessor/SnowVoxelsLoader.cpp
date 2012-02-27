/*************************************************************************
*  Copyright (C) 2008 by Janek Kozicki <cosurgi@berlios.de>              *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

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
#include<yade/pkg-dem/TriaxialStateRecorder.hpp>

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

#include<yade/pkg-snow/ElawSnowLayersDeformation.hpp>

#include<yade/pkg-dem/Shop.hpp>

#include<boost/filesystem/convenience.hpp>
#include<boost/lexical_cast.hpp>
#include<boost/numeric/conversion/bounds.hpp>
#include<boost/limits.hpp>


#include"SnowVoxelsLoader.hpp"
#include<set>
#include<yade/pkg-snow/Ef2_BssSnowGrain_BssSnowGrain_makeIstSnowLayersContact.hpp>
#include<yade/pkg-snow/Ef2_InteractingBox_BssSnowGrain_makeIstSnowLayersContact.hpp>
//#include<yade/pkg-snow/Ef2_BssSnowGrain_BssSnowGrain_makeSpheresContactGeometry.hpp>
//#include<yade/pkg-snow/Ef2_InteractingBox_BssSnowGrain_makeSpheresContactGeometry.hpp>

/* this is a little hack; all those files are not picked up by the source
scanner since they have YADE_REQUIRE_FEATURE(some nonsense) inside.
Sorry for that. */
#include"Voxel/Config.cpp"
#include"Voxel/DataSurface.cpp"
#include"Voxel/DataVoxel.cpp"
#include"Voxel/GrainSurface.cpp"
#include"Voxel/SafeVectors3.cpp"
#include"Voxel/VoxelEnvelope.cpp"

YADE_REQUIRE_FEATURE(geometricalmodel);
YADE_PLUGIN((SnowVoxelsLoader));

SnowVoxelsLoader::SnowVoxelsLoader() : FileGenerator()
{
	voxel_binary_data_file = "/home/janek/32-Snow-white/20-Programy/31-SNOW-read-data/RESULT.bz2";
//	voxel_txt_dir = "";
//	voxel_caxis_file = "";
//	voxel_colors_file = "";
	grain_binary_data_file = "grain_binary--NOT-USED-NOW.bz2";
	
	sphereYoungModulus	= 15000000.0;
	spherePoissonRatio	= 0.5;
	sphereFrictionDeg	= 18.0;
	boxYoungModulus		= 15000000.0;
	boxPoissonRatio		= 0.2;
	boxFrictionDeg		= 0.f;
	density			= 1100;
	
	dampingForce		= 0.2;
	dampingMomentum		= 0.2;
	defaultDt		= 0.00001;
	timeStepUpdateInterval	= 50;
	wallStiffnessUpdateInterval = 10;
	radiusControlInterval = 10;

	spheresColor		= Vector3r(0.8,0.3,0.3);
	use_grain_shear_creep	= true;
	use_grain_twist_creep	= true;

// a pixel is 20.4 microns (2.04 � 10-5 meters)
// the sample was 10.4mm high
	one_voxel_in_meters_is	= 2.04e-5;
	layer_distance_voxels   = 6.0;
	angle_increment_radians = 0.3;
	skip_small_grains       = 25000;

	WallStressRecordFile	= "./SnowWallStresses";

	normalCohesion		= 50000000;
	shearCohesion		= 50000000;
	setCohesionOnNewContacts= true;
	
	sigma_iso = 50000;
	creep_viscosity = 4000000;
	
	thickness 		= 0.0003;
	strainRate = 10;
	StabilityCriterion = 0.01;
	autoCompressionActivation = false;
	internalCompaction	=false;
	maxMultiplier = 1.01;
	
	gravity			= Vector3r(0,0,-9.81);
	
	recordIntervalIter	= 20;

	use_gravity_engine = false;
	enable_layers_creep = true;
	layers_creep_viscosity = 100000;
	m_grains.clear();
}


SnowVoxelsLoader::~SnowVoxelsLoader ()
{
}

bool SnowVoxelsLoader::load_voxels()
{
//	if(grain_binary_data_file !="" && boost::filesystem::exists(grain_binary_data_file))
//	{
//		std::cerr << "no need to load voxels - grain binary file exists\n";
//		std::cerr << "loading " << grain_binary_data_file << " ...";
//			boost::iostreams::filtering_istream ifs;
//			ifs.push(boost::iostreams::bzip2_decompressor());
//			ifs.push(boost::iostreams::file_source(grain_binary_data_file));
//		//std::ifstream ifs(m_config.load_file());
//		boost::archive::binary_iarchive ia(ifs);
//		ia >> m_grains;
//		return true;
//	}

	if(voxel_binary_data_file != "" && boost::filesystem::exists(voxel_binary_data_file))
	{
		std::cerr << "loading " << voxel_binary_data_file << " ...";
			boost::iostreams::filtering_istream ifs;
			ifs.push(boost::iostreams::bzip2_decompressor());
			ifs.push(boost::iostreams::file_source(voxel_binary_data_file));
		//std::ifstream ifs(m_config.load_file());
		boost::archive::binary_iarchive ia(ifs);
		ia >> m_voxel;
		std::cerr << " finished\n";
		return true;
	}
	else
	{
		message="cannot load input file check if voxel_binary_data_file exists. Or check if grain_binary_data_file does not exist - because if it exists it is loaded first, skipping the layers generation.";
		return false;
	}
	return false;
}

int voxel_id_count(int id,const T_DATA& dat)
{	
	int res=0;	
	BOOST_FOREACH(const std::vector<std::vector<unsigned char> >& a,dat)
		BOOST_FOREACH(const std::vector<unsigned char>& b,a)
			BOOST_FOREACH(unsigned char c,b)
				if(c==id)
					++res;
	return res;
}

bool SnowVoxelsLoader::generate()
{
	if(!load_voxels())
		return false;

	rootBody = shared_ptr<Scene>(new Scene);
	createActors(rootBody);
	positionRootBody(rootBody);

	
	if(m_grains.size() == 0)
	{
		int skip_total(0);
		const T_DATA& dat(m_voxel.m_data.get_data_voxel_const_ref().get_a_voxel_const_ref());
		std::set<unsigned char> done;done.clear();done.insert(0);
		BOOST_FOREACH(const std::vector<std::vector<unsigned char> >& a,dat)
			BOOST_FOREACH(const std::vector<unsigned char>& b,a)
				BOOST_FOREACH(unsigned char c,b)
				{
					if(done.find(c)==done.end())
					{
						done.insert(c);
						if(voxel_id_count(c,dat) > skip_small_grains )
						{
							boost::shared_ptr<BshSnowGrain> grain(new BshSnowGrain(dat,m_voxel.m_data.get_axes_const_ref()[c],c,m_voxel.m_data.get_colors_const_ref()[c],one_voxel_in_meters_is,layer_distance_voxels,angle_increment_radians));
							m_grains.push_back(grain);
						}
						else
						// skip too small grains, because they require extremely small timestep, 
						// and calculations will take forever.
						// (they bounce a lot with bigger timestep)
						{
							std::cerr << "\n======= skipped grain id: " << ((int)(c)) << "\n";
							++skip_total;
						}
					}
				};
		
		std::cerr << "\n======= total skipped grains: " << ((int)(skip_total)) << "\n";

//		std::cerr << "saving "<< grain_binary_data_file << " ...";
//		boost::iostreams::filtering_ostream ofs;
//		ofs.push(boost::iostreams::bzip2_compressor());
//		ofs.push(boost::iostreams::file_sink(grain_binary_data_file));
//		boost::archive::binary_oarchive oa(ofs);
////#if BOOST_VERSION >= 103500
//		oa << m_grains;
//#else
//		const std::vector<boost::shared_ptr<BshSnowGrain> > tmp(m_grains);
//		oa << tmp;
//#endif
		std::cerr << " finished\n";
	}
	
	Vector3r upperCorner(-Mathr::MAX_REAL,-Mathr::MAX_REAL,-Mathr::MAX_REAL);
	Vector3r lowerCorner( Mathr::MAX_REAL, Mathr::MAX_REAL, Mathr::MAX_REAL);
	
	shared_ptr<Body> body;
	BOOST_FOREACH(boost::shared_ptr<BshSnowGrain> grain,m_grains)
	{
		//lowerCorner, upperCorner
		BOOST_FOREACH(std::vector<Vector3r>& g,grain->slices)
			BOOST_FOREACH(Vector3r& v,g)
			{
				upperCorner = componentMaxVector(upperCorner,v+grain->center);
				lowerCorner = componentMinVector(lowerCorner,v+grain->center);
			}
	}
	Real sx = upperCorner[0] - lowerCorner[0];
	Real sy = upperCorner[1] - lowerCorner[1];
// make box a little bigger - the four walls are away 10% of size x (sx) and size y (sy)
	upperCorner[0]+=sx*0.1;
	lowerCorner[0]-=sx*0.1;
	upperCorner[1]+=sy*0.1;
	lowerCorner[1]-=sy*0.1;
//	thickness = sx*0.05;
	// Z touches top and bottom, exactly.
	// but X,Y give some free space
	{ // Now make 6 walls for the TriaxialCompressionEngine
	// bottom box
	 	Vector3r center		= Vector3r(
	 						(lowerCorner[0]+upperCorner[0])/2,
	 						lowerCorner[1]-thickness/2.0,
	 						(lowerCorner[2]+upperCorner[2])/2);
	 	Vector3r halfSize	= Vector3r(
	 						1.5*fabs(lowerCorner[0]-upperCorner[0])/2+thickness,
							thickness/2.0,
	 						1.5*fabs(lowerCorner[2]-upperCorner[2])/2+thickness);
	
		create_box(body,center,halfSize,true);
			rootBody->bodies->insert(body);
			triaxialcompressionEngine->wall_bottom_id = body->getId();
	
	// top box
	 	center			= Vector3r(
	 						(lowerCorner[0]+upperCorner[0])/2,
	 						upperCorner[1]+thickness/2.0,
	 						(lowerCorner[2]+upperCorner[2])/2);
	 	halfSize		= Vector3r(
	 						1.5*fabs(lowerCorner[0]-upperCorner[0])/2+thickness,
	 						thickness/2.0,
	 						1.5*fabs(lowerCorner[2]-upperCorner[2])/2+thickness);
	
		create_box(body,center,halfSize,true);
			rootBody->bodies->insert(body);
			triaxialcompressionEngine->wall_top_id = body->getId();
	// box 1
	
	 	center			= Vector3r(
	 						lowerCorner[0]-thickness/2.0,
	 						(lowerCorner[1]+upperCorner[1])/2,
	 						(lowerCorner[2]+upperCorner[2])/2);
		halfSize		= Vector3r(
							thickness/2.0,
	 						1.5*fabs(lowerCorner[1]-upperCorner[1])/2+thickness,
	 						1.5*fabs(lowerCorner[2]-upperCorner[2])/2+thickness);
		create_box(body,center,halfSize,true);
			rootBody->bodies->insert(body);
			triaxialcompressionEngine->wall_left_id = body->getId();
	// box 2
	 	center			= Vector3r(
	 						upperCorner[0]+thickness/2.0,
	 						(lowerCorner[1]+upperCorner[1])/2,
							(lowerCorner[2]+upperCorner[2])/2);
	 	halfSize		= Vector3r(
	 						thickness/2.0,
	 						1.5*fabs(lowerCorner[1]-upperCorner[1])/2+thickness,
	 						1.5*fabs(lowerCorner[2]-upperCorner[2])/2+thickness);
	 	
		create_box(body,center,halfSize,true);
			rootBody->bodies->insert(body);
			triaxialcompressionEngine->wall_right_id = body->getId();
	// box 3
	 	center			= Vector3r(
	 						(lowerCorner[0]+upperCorner[0])/2,
	 						(lowerCorner[1]+upperCorner[1])/2,
	 						lowerCorner[2]-thickness/2.0);
	 	halfSize		= Vector3r(
	 						1.5*fabs(lowerCorner[0]-upperCorner[0])/2+thickness,
	 						1.5*fabs(lowerCorner[1]-upperCorner[1])/2+thickness,
	 						thickness/2.0);
		create_box(body,center,halfSize,true);
			rootBody->bodies->insert(body);
			triaxialcompressionEngine->wall_back_id = body->getId();
	
	// box 4
	 	center			= Vector3r(
	 						(lowerCorner[0]+upperCorner[0])/2,
	 						(lowerCorner[1]+upperCorner[1])/2,
	 						upperCorner[2]+thickness/2.0);
	 	halfSize		= Vector3r(
	 						1.5*fabs(lowerCorner[0]-upperCorner[0])/2+thickness,
	 						1.5*fabs(lowerCorner[1]-upperCorner[1])/2+thickness,
	 						thickness/2.0);
		create_box(body,center,halfSize,true);
			rootBody->bodies->insert(body);
			triaxialcompressionEngine->wall_front_id = body->getId();
			 
	}
	BOOST_FOREACH(boost::shared_ptr<BshSnowGrain> grain,m_grains)
	{
		create_grain(body,grain->center,true,grain);
		rootBody->bodies->insert(body);
	}
	
	return true;
}

void SnowVoxelsLoader::createActors(shared_ptr<Scene>& rootBody)
{
	shared_ptr<InteractionGeometryDispatcher> interactionGeometryDispatcher(new InteractionGeometryDispatcher);

	shared_ptr<InteractionGeometryFunctor> s1(new Ef2_BssSnowGrain_BssSnowGrain_makeIstSnowLayersContact);
	shared_ptr<InteractionGeometryFunctor> s2(new Ef2_InteractingBox_BssSnowGrain_makeIstSnowLayersContact);
	//shared_ptr<InteractionGeometryFunctor> s1(new Ef2_BssSnowGrain_BssSnowGrain_makeSpheresContactGeometry);
	//shared_ptr<InteractionGeometryFunctor> s2(new Ef2_InteractingBox_BssSnowGrain_makeSpheresContactGeometry);
	//shared_ptr<InteractionGeometryFunctor> s1(new Ig2_Sphere_Sphere_ScGeom);
	//shared_ptr<InteractionGeometryFunctor> s2(new Ig2_Box_Sphere_ScGeom);
	
	interactionGeometryDispatcher->add(s1);
	interactionGeometryDispatcher->add(s2);

	shared_ptr<Ip2_2xCohFrictMat_CohFrictPhys> cohesiveFrictionalRelationships = shared_ptr<Ip2_2xCohFrictMat_CohFrictPhys> (new Ip2_2xCohFrictMat_CohFrictPhys);
	cohesiveFrictionalRelationships->shearCohesion = shearCohesion;
	cohesiveFrictionalRelationships->normalCohesion = normalCohesion;
	cohesiveFrictionalRelationships->setCohesionOnNewContacts = setCohesionOnNewContacts;
	cohesiveFrictionalRelationships->setCohesionNow = true;
	shared_ptr<InteractionPhysicsDispatcher> interactionPhysicsDispatcher(new InteractionPhysicsDispatcher);
	interactionPhysicsDispatcher->add(cohesiveFrictionalRelationships);
		
	shared_ptr<BoundDispatcher> boundDispatcher	= shared_ptr<BoundDispatcher>(new BoundDispatcher);
	boundDispatcher->add("Ef2_BssSnowGrain_AABB_makeAABB");
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

	shared_ptr<GlobalStiffnessTimeStepper> globalStiffnessTimeStepper(new GlobalStiffnessTimeStepper);
	globalStiffnessTimeStepper->sdecGroupMask = 2;
	globalStiffnessTimeStepper->timeStepUpdateInterval = timeStepUpdateInterval;
	globalStiffnessTimeStepper->defaultDt = defaultDt;
	globalStiffnessTimeStepper->timestepSafetyCoefficient = 0.1;
	
	shared_ptr<CohesiveFrictionalContactLaw> cohesiveFrictionalContactLaw(new CohesiveFrictionalContactLaw);
	cohesiveFrictionalContactLaw->sdecGroupMask = 2;
	cohesiveFrictionalContactLaw->shear_creep = use_grain_shear_creep;
	cohesiveFrictionalContactLaw->twist_creep = use_grain_twist_creep;
	cohesiveFrictionalContactLaw->creep_viscosity = creep_viscosity;
		
	
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
	triaxialStateRecorder = shared_ptr<TriaxialStateRecorder>(new TriaxialStateRecorder);
	triaxialStateRecorder-> outputFile 	= WallStressRecordFile;
	triaxialStateRecorder-> interval 		= recordIntervalIter;
	//triaxialStateRecorder-> thickness 		= thickness;
	
	shared_ptr<ElawSnowLayersDeformation>elawSnowLayersDeformation(new ElawSnowLayersDeformation);
	elawSnowLayersDeformation->creep_viscosity = layers_creep_viscosity;
	elawSnowLayersDeformation->sdecGroupMask = 2;

	// moving walls to regulate the stress applied
	//cerr << "triaxialstressController = shared_ptr<TriaxialStressController> (new TriaxialStressController);" << std::endl;
	triaxialstressController = YADE_PTR_CAST<TriaxialStressController>(triaxialcompressionEngine);
	triaxialstressController-> stiffnessUpdateInterval = wallStiffnessUpdateInterval;// = recordIntervalIter
	triaxialstressController-> sigma_iso = sigma_iso;
	triaxialstressController-> max_vel = 1;
	triaxialstressController-> thickness = thickness;
	triaxialstressController->wall_bottom_activated = false;
	triaxialstressController->wall_top_activated = false;
	triaxialstressController->wall_left_activated = false;
	triaxialstressController->wall_right_activated = false;
	triaxialstressController->wall_front_activated = true;
	triaxialstressController->wall_back_activated = true;
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
	if(use_gravity_engine)
		rootBody->engines.push_back(gravityCondition);
	rootBody->engines.push_back(actionDampingDispatcher);
	if(enable_layers_creep)
		rootBody->engines.push_back(elawSnowLayersDeformation);
	rootBody->engines.push_back(applyActionDispatcher);
	rootBody->engines.push_back(positionIntegrator);
	//if(!rotationBlocked)
		rootBody->engines.push_back(orientationIntegrator);
	//rootBody->engines.push_back(triaxialstressController);
	
		
	//rootBody->engines.push_back(averagePositionRecorder);
	//rootBody->engines.push_back(velocityRecorder);
	//rootBody->engines.push_back(forcerec);
	
	rootBody->initializers.clear();
	rootBody->initializers.push_back(boundDispatcher);
	
}

void SnowVoxelsLoader::positionRootBody(shared_ptr<Scene>& rootBody)
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

void SnowVoxelsLoader::create_grain(shared_ptr<Body>& body, Vector3r position, bool dynamic , boost::shared_ptr<BshSnowGrain> grain)
{
	body = shared_ptr<Body>(new Body(body_id_t(0),2));
	shared_ptr<CohFrictMat> physics(new CohFrictMat);
	shared_ptr<Aabb> aabb(new Aabb);
	shared_ptr<BshSnowGrain> gSnowGrain(grain);
	

	shared_ptr<BssSnowGrain> iSphere(new BssSnowGrain(gSnowGrain.get(),one_voxel_in_meters_is)); 
	Real radius = iSphere->radius;
	//shared_ptr<Sphere> iSphere(new Sphere);
	//Real radius = (grain->start-grain->end).Length()*0.5;
	
	Quaternionr q;//(Mathr::SymmetricRandom(),Mathr::SymmetricRandom(),Mathr::SymmetricRandom(),Mathr::SymmetricRandom());
	q.FromAxisAngle( Vector3r(0,0,1),0);
	q.Normalize();
	
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

	aabb->diffuseColor		= Vector3r(0,1,0);


	gSnowGrain->diffuseColor	= grain->color;
	gSnowGrain->wire		= false;
	gSnowGrain->shadowCaster	= true;
	
	//iSphere->radius			= radius; // already calculated
	iSphere->diffuseColor		= grain->color;

	body->shape	= iSphere;
	body->geometricalModel		= gSnowGrain;
	body->bound		= aabb;
	body->physicalParameters	= physics;
}


void SnowVoxelsLoader::create_box(shared_ptr<Body>& body, Vector3r position, Vector3r extents, bool wire)
{
	body = shared_ptr<Body>(new Body(body_id_t(0),2));
	shared_ptr<CohFrictMat> physics(new CohFrictMat);
	shared_ptr<Aabb> aabb(new Aabb);
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

	
	iBox->extents			= extents;
	iBox->diffuseColor		= Vector3r(0.5,0.5,0.5);

	body->bound		= aabb;
	body->shape	= iBox;
	#ifdef YADE_GEOMETRICALMODEL
		shared_ptr<BoxModel> gBox(new BoxModel);
		gBox->extents			= extents;
		gBox->diffuseColor		= Vector3r(0.5,0.5,0.5);
		gBox->wire			= wire;
		gBox->shadowCaster		= false;
		body->geometricalModel		= gBox;
	#endif
	body->physicalParameters	= physics;
}


YADE_REQUIRE_FEATURE(PHYSPAR);

