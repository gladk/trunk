/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include<iostream>
#include<fstream>
#include"LatticeExampleCTData.hpp"
#include<yade/core/Scene.hpp>
#include<yade/pkg-common/BoundDispatcher.hpp>
#include<yade/pkg-common/GeometricalModelMetaEngine.hpp>
#include<yade/pkg-common/Aabb.hpp>
#include<yade/pkg-common/SphereModel.hpp>
#include<yade/pkg-common/DisplacementEngine.hpp>
#include<yade/pkg-lattice/LatticeLaw.hpp>
#include<yade/pkg-lattice/LatticeSetParameters.hpp>
#include<yade/pkg-lattice/LatticeInteractingGeometry.hpp>
#include<yade/pkg-lattice/LatticeNodeParameters.hpp>
#include<yade/pkg-lattice/LatticeBeamParameters.hpp>
#include<yade/pkg-lattice/LineSegment.hpp>
#include<yade/pkg-lattice/LatticeBeamAngularSpring.hpp>

// Delaunay
#ifndef MINIWM3
	#include <Wm3Delaunay3.h>
	#include <Wm3Delaunay2.h>
	#include <Wm3Query.h>
#endif

using namespace boost;
using namespace std;

YADE_REQUIRE_FEATURE(geometricalmodel);

LatticeExampleCTData::LatticeExampleCTData() : FileGenerator()
{
	nodeGroupMask                   = 1;
	beamGroupMask                   = 2;

	// mesh generation
	speciemen_size_in_meters        = Vector3r(0.1,0.1,0.1);
	cell_size_in_meters             = 0.02;
	//disorder_in_cellsize_unit       = Vector3r(0.6,0.6,0.6);
	disorder_in_cellsize_unit       = Vector3r(0.0,0.0,0.0);
	//max_rod_length_in_cellsize_unit = 1.5;
	max_rod_length_in_cellsize_unit = 1.01;

	// MaterialParameters of cement matrix
	crit_tensile_strain             = 1.0;       // E_min
	crit_compr_strain               = 0.5;       // E_max

	longitudal_stiffness_no_unit    = 1.0;        // k_l
	bending_stiffness_no_unit       = 0.6;        // k_b
	torsional_stiffness_no_unit     = 0.6;        // k_t

	// conditions - displacement of top & bottom
	smooth_edges_displacement       = false;
	y_top_bottom_displ_in_meters    = 0.00001;

	ct_data_file_name               = "usg_rainier_8_2_1_50_20LH_128x128x128_binary.raw";
}


LatticeExampleCTData::~LatticeExampleCTData()
{

}

bool LatticeExampleCTData::generate()
{
	vector<vector<vector<unsigned char> > > ct_data;
	ct_data.resize(128);
	BOOST_FOREACH(vector<vector<unsigned char> >& b,ct_data)
	{
		b.resize(128);
		BOOST_FOREACH(vector<unsigned char>& c,b)
			c.resize(128);
	}
	ifstream file(ct_data_file_name.c_str(),ios::in | ios::binary);
	if(!file)
	{
		message=ct_data_file_name + " does not exist! Cannot load CT data.";
		return false;
	}
	for(int i=0;i<128;i++)
		for(int j=0;j<128;j++)
		{
			for(int k=0;k<128;k++)
			{
				file.read((char *)(&ct_data[i][j][k]),1);
				if(ct_data[i][j][k] > 10)
					std::cout << "#";
				else
					std::cout << ".";
//				std::cerr << (int)(ct_data[i][j][k]) << "\n";
			}
		std::cout << "\n";
	}
//	return false;


	rootBody = shared_ptr<Scene>(new Scene);
	make_simulation_loop(rootBody);    // make the simulation loop
	positionRootBody(rootBody); // set global coordinate system, etc.
	

	shared_ptr<Body> body;

// the number of nodes to generate in each direction
	Vector3r nbNodes = speciemen_size_in_meters / cell_size_in_meters + Vector3r(1,1,1);
	std::vector<Vector3r> vert3; vert3.clear(); // Delaunay 3D

// create nodes
	setStatus("creating nodes...");
        unsigned int totalNodesCount = 0;
	unsigned int totalBeamsCount = 0;
	float all = nbNodes[0]*nbNodes[1]*nbNodes[2];
	float current = 0.0;

// loop on each node in a grid, to create it
	for( int j=0 ; j<=nbNodes[1] ; j++ )
		for( int i=0 ; i<=nbNodes[0] ; i++ )
			for( int k=0 ; k<=nbNodes[2] ; k++)
			{
				shared_ptr<Body> node;

				// calculate node position - orthogonal grid, see LatticeExample for tetrahedral grid
				// WARNING - orthogonal grid may produce unrealistic results!
				Vector3r position = Vector3r(i,j,k) * cell_size_in_meters;
				position += Vector3r( 	  Mathr::SymmetricRandom()*disorder_in_cellsize_unit[0]
					  		, Mathr::SymmetricRandom()*disorder_in_cellsize_unit[1]
							, Mathr::SymmetricRandom()*disorder_in_cellsize_unit[2]) * cell_size_in_meters;

				int I = (int)(128.0*position[0]/speciemen_size_in_meters[0]);
				int J = (int)(128.0*position[1]/speciemen_size_in_meters[1]);
				int K = (int)(128.0*position[2]/speciemen_size_in_meters[2]);
				I = std::max(0,I); I=std::min(I,127);
				J = std::max(0,J); J=std::min(J,127);
				K = std::max(0,K); K=std::min(K,127);
				if(ct_data[I][J][K] == 0)
				{ // check if node position is not the air

					if(createNodeXYZ(node,position))
					{
						rootBody->bodies->insert(node);
						++totalNodesCount;
						vert3.push_back(node->physicalParameters->se3.position); // Delaunay 3D
					}
				}
				
				setProgress(current++/all);
				if(shouldTerminate()) return false;
			}


// now connect the nodes with rods
	setStatus("creating rods...");

#ifdef MINIWM3
// if Delaunay is not avauilable produce an error message
	message = "ERROR:\n\nYade was compiled without full installation of wildmagic-dev (Wm3 foundation library), can't use Delaunay. \n"
	"Please compile yade with full wm3 library installed. Don't forget CPPPATH=/usr/include/wm3 scons parameter.";
	return false;
#else
// create beams, Delaunay
	std::set< std::pair<int,int> > pairs;
	int I,J;
	setProgress(0); setStatus("Delaunay 3d...");
	Delaunay3<Real> del3(vert3.size(),&(vert3[0]),cell_size_in_meters/50,false, Query::QT_INTEGER); // Delaunay
	
	int del3_i = 0;
	int del3_ind[4];
	while(del3.GetIndexSet(del3_i++,del3_ind))
	{
		if(shouldTerminate()) return false;
		setProgress((float)del3_i/(float)del3.GetSimplexQuantity());
		shared_ptr<Body> beam;

		for(int delauney3d_index=0 ; delauney3d_index<4 ; ++delauney3d_index)
		{
			I=del3_ind[delauney3d_index];J=del3_ind[(delauney3d_index+1)%4];
			if(J>I) std::swap(I,J);
			if(pairs.insert(std::make_pair(I,J)).second)
			{
				createBeam(beam,I,J);

				//
				// pos1 and pos2 are nodes' positions for this beam, might be useful if beam 
				// parameters depend on position in the CT data
				//
				Vector3r pos1 = (*(rootBody->bodies))[I]->physicalParameters->se3.position;
				Vector3r pos2 = (*(rootBody->bodies))[J]->physicalParameters->se3.position;
		
				////////////////////////////////////////////////////////////////////////
				// the beam parameters.
				////////////////////////////////////////////////////////////////////////
				Real E_min = crit_tensile_strain;
				Real E_max = crit_compr_strain;
				Real k_l   = longitudal_stiffness_no_unit;
				Real k_b   = bending_stiffness_no_unit;
				Real k_t   = torsional_stiffness_no_unit;
				////////////////////////////////////////////////////////////////////////
				// the beam parameters END.
				////////////////////////////////////////////////////////////////////////

				// insert rod, only if it is short enough
				if(calcBeamPositionOrientationLength(beam,E_min,E_max,k_l,k_b,k_t) < max_rod_length_in_cellsize_unit*cell_size_in_meters)
				{
					rootBody->bodies->insert(beam);
					++totalBeamsCount;
				}
			}
		}
	}
#endif

	create_angular_springs(totalNodesCount,rootBody);
	if(shouldTerminate()) return false;

// calculate the displacement regions, and assign displacements
	Real disp = y_top_bottom_displ_in_meters;
	Vector3r dir(0,1,0);
	Vector3r min1(-1,-1,-1);
	Vector3r max1(speciemen_size_in_meters[0]+1 , cell_size_in_meters*0.8 , speciemen_size_in_meters[2]+1 );
	Vector3r min2(-1, speciemen_size_in_meters[1]-cell_size_in_meters*1.1 , -1);
	Vector3r max2(speciemen_size_in_meters*2);

	imposeTranslation(rootBody,min1,max1,-1.0 * dir, disp);
	imposeTranslation(rootBody,min2,max2,       dir, disp);

	if(shouldTerminate()) return false;

        message="Number of nodes created: "  + lexical_cast<string>(totalNodesCount) +
					" (" + lexical_cast<string>(nbNodes[0]) + " x "
					     + lexical_cast<string>(nbNodes[1]) + " x "
					     + lexical_cast<string>(nbNodes[2]) + ")\n\n"
	     + "Number of beams: " + lexical_cast<string>(totalBeamsCount) + "\n\n"
	     + "NOTE: sometimes it can look better when 'drawWireFrame' is enabled in Display tab.\n";

	message = "WARNING - orthogonal grid may produce unrealistic results! See LatticeExample.cpp for tetrahedral grid.\n\n" + message;

        cerr << "finished.. saving\n" << message << "\n";

	return true;
}

void LatticeExampleCTData::create_angular_springs(int totalNodesCount,shared_ptr<Scene> rootBody)
{
	// check what nodes are connected with which beams, this is necessary for creating the angular/torsional springs
        std::vector< std::vector< unsigned int > > connections;

	BodyContainer::iterator bi    = rootBody->bodies->begin();
	BodyContainer::iterator biEnd = rootBody->bodies->end();
	connections.resize(totalNodesCount);
	for(  ; bi!=biEnd ; ++bi )  // loop over all beams
	{
		if(shouldTerminate()) return;

		Body* body = (*bi).get();
		if( ! ( body->getGroupMask() & beamGroupMask ) )
			continue; // skip non-beams
		
		LatticeBeamParameters* beam = static_cast<LatticeBeamParameters*>(body->physicalParameters.get() );
		connections[beam->id1].push_back(body->getId());
		connections[beam->id2].push_back(body->getId());
	}

	// create angular springs between beams
	setStatus("angular springs...");
	bi    = rootBody->bodies->begin();
	biEnd = rootBody->bodies->end();
	float all_bodies = rootBody->bodies->size();
	int current = 0;
	for(  ; bi!=biEnd ; ++bi )  // loop over all beams
	{
		if(shouldTerminate()) return;

		if( ++current % 100 == 0 )
		{
			//cerr << "angular springs: " << current << " , " << ((static_cast<float>(current)/all_bodies)*100.0) << " %\n";
			setProgress(((float)(current)/all_bodies));
		}

		Body* body = (*bi).get();
		if( ! ( body->getGroupMask() & beamGroupMask ) )
			continue; // skip non-beams

		calcBeamAngles(body,rootBody->bodies.get(),rootBody->interactions.get(),connections);
	}
}

bool LatticeExampleCTData::createNodeXYZ(shared_ptr<Body>& body, Vector3r position)
{
	body = shared_ptr<Body>(new Body(body_id_t(0),nodeGroupMask));
	shared_ptr<LatticeNodeParameters> physics(new LatticeNodeParameters);
	shared_ptr<SphereModel> gSphere(new SphereModel);
	
	Quaternionr q;
	q.FromAxisAngle( Vector3r(Mathr::UnitRandom(),Mathr::UnitRandom(),Mathr::UnitRandom()) , Mathr::UnitRandom()*Mathr::PI );
	
	Real radius 			= cell_size_in_meters*0.05;
	
	body->isDynamic			= true;
	
	physics->se3			= Se3r(position,q);

	gSphere->radius			= radius;
	gSphere->diffuseColor		= Vector3r(0.8,0.8,0.8);
	gSphere->wire			= false;
	gSphere->shadowCaster		= false;
	
        body->geometricalModel          = gSphere;
        body->physicalParameters        = physics;
 
	if( 	   position[0] >= speciemen_size_in_meters[0] 
		|| position[1] >= speciemen_size_in_meters[1]
		|| position[2] >= speciemen_size_in_meters[2] )
		return false;
	if( 	   position[0] < 0 
		|| position[1] < 0
		|| position[2] < 0 )
		return false;

        return true;
}


void LatticeExampleCTData::createBeam(shared_ptr<Body>& body, int i, int j)
{
	body = shared_ptr<Body>(new Body(body_id_t(0),beamGroupMask));
	shared_ptr<LatticeBeamParameters> physics(new LatticeBeamParameters);
	shared_ptr<LineSegment> gBeam(new LineSegment);
	
	Real length 			= 1.0; // unspecified for now, calcBeamPositionOrientationLength will calculate it
	
	body->isDynamic			= true;
	
	physics->id1 			= i;
	physics->id2 			= j;

	gBeam->length			= length;
	gBeam->diffuseColor		= Vector3r(0.6,0.6,0.6);
	gBeam->wire			= false;
	gBeam->shadowCaster		= false;
	
	body->geometricalModel		= gBeam;
	body->physicalParameters	= physics;
}


Real LatticeExampleCTData::calcBeamPositionOrientationLength(
	shared_ptr<Body>& body,
	Real E_min,
	Real E_max,
	Real k_l,
	Real k_b,
	Real k_t)
{
	LatticeBeamParameters* beam = dynamic_cast<LatticeBeamParameters*>(body->physicalParameters.get());
	LineSegment* ls = dynamic_cast<LineSegment*>(body->geometricalModel.get());
	if(beam != 0 && ls != 0)
	{
		shared_ptr<Body>& bodyA = (*(rootBody->bodies))[beam->id1];
		shared_ptr<Body>& bodyB = (*(rootBody->bodies))[beam->id2];
		Se3r& se3A 		= bodyA->physicalParameters->se3;
		Se3r& se3B 		= bodyB->physicalParameters->se3;
		
		Se3r se3Beam;
		se3Beam.position 	= (se3A.position + se3B.position)*0.5;
		Vector3r dist 		= se3A.position - se3B.position;
		
		Real length 		= dist.Normalize();
		beam->direction 	= dist;
		beam->length            = length;
		beam->initialLength     = length;
		ls->length = length;

		beam->criticalTensileStrain     = E_min;
		beam->criticalCompressiveStrain = E_max;
		beam->longitudalStiffness       = k_l;
		beam->bendingStiffness          = k_b;
		beam->torsionalStiffness        = k_t;
		
		se3Beam.orientation.Align( Vector3r::UnitX() , dist );
		beam->se3 		= se3Beam;
		beam->se3Displacement.position 	= Vector3r(0.0,0.0,0.0);
		beam->se3Displacement.orientation.Align(dist,dist);

		beam->otherDirection	= beam->se3.orientation*Vector3r::UnitY(); // any unit vector that is orthogonal to direction.

		return length;
	}
	else
	{
		std::cerr << body->getId() << " is not a beam, skipping\n";
		return 10000000;
	}
}

void LatticeExampleCTData::calcAxisAngle(LatticeBeamParameters* beam1, BodyContainer* bodies, int otherId, InteractionContainer* ints, int thisId)
{ 
	if( ! ints->find(body_id_t(otherId),body_id_t(thisId)) && otherId != thisId )
	{
		LatticeBeamParameters* 	beam2 		= static_cast<LatticeBeamParameters*>( ((*(bodies))[ otherId ])->physicalParameters.get() );
		Real 			angle;
		
		angle = unitVectorsAngle(beam1->direction,beam2->direction);

                shared_ptr<Interaction>                 interaction(new Interaction( body_id_t(thisId) , body_id_t(otherId) ));
                shared_ptr<LatticeBeamAngularSpring>    angularSpring(new LatticeBeamAngularSpring);
		
		angularSpring->initialPlaneAngle 	= angle;
		angularSpring->planeSwap180		= false;
		angularSpring->lastCrossProduct		= 1.0*(beam1->direction.Cross(beam2->direction));
		angularSpring->lastCrossProduct.Normalize();
		
		angularSpring->initialOffPlaneAngle1	= unitVectorsAngle(beam1->otherDirection,angularSpring->lastCrossProduct);
		angularSpring->initialOffPlaneAngle2	= unitVectorsAngle(beam2->otherDirection,angularSpring->lastCrossProduct);

	//	std::cerr << thisId << ", " << otherId << ", (" << beam1->otherDirection << ") (" << beam2->otherDirection << ") (" << angularSpring->lastCrossProduct << ")\n";
				
		Quaternionr	aligner1,aligner2;
		aligner1.FromAxisAngle(beam1->direction , angularSpring->initialOffPlaneAngle1);
		aligner2.FromAxisAngle(beam2->direction , angularSpring->initialOffPlaneAngle2);
				
		Vector3r	dir1			= aligner1 * angularSpring->lastCrossProduct;
		Vector3r	dir2			= aligner2 * angularSpring->lastCrossProduct;

		// insignificant error is possible here
		if( dir1.Dot(beam1->otherDirection) < 0.999999 )
			angularSpring->initialOffPlaneAngle1   *= -1.0;//, angularSpring->offPlaneSwap1 = true;
		if( dir2.Dot(beam2->otherDirection) < 0.999999 )
			angularSpring->initialOffPlaneAngle2   *= -1.0;//, angularSpring->offPlaneSwap2 = true;
	
		interaction->interactionPhysics 	= angularSpring;
		ints->insert(interaction);
	}
}

void LatticeExampleCTData::calcBeamAngles(
	Body* body, 
	BodyContainer* bodies, 
	InteractionContainer* ints,
	std::vector< std::vector< unsigned int > >& connections)
{
	LatticeBeamParameters* beam 	= static_cast<LatticeBeamParameters*>(body->physicalParameters.get());

	std::vector<unsigned int>::iterator i   = connections[beam->id1].begin();
	std::vector<unsigned int>::iterator end = connections[beam->id1].end();
	
	for( ; i != end ; ++i )
		calcAxisAngle(beam,bodies,*i,ints,body->getId());
	
	i   = connections[beam->id2].begin();
	end = connections[beam->id2].end();
	for( ; i != end ; ++i )
                calcAxisAngle(beam,bodies,*i,ints,body->getId());
}


void LatticeExampleCTData::make_simulation_loop(shared_ptr<Scene>& )
{
	shared_ptr<BoundDispatcher> boundDispatcher   = shared_ptr<BoundDispatcher>(new BoundDispatcher);

	shared_ptr<GeometricalModelMetaEngine> geometricalModelDispatcher       = shared_ptr<GeometricalModelMetaEngine>(new GeometricalModelMetaEngine);
	geometricalModelDispatcher->add("LatticeSet2LatticeBeams");

	shared_ptr<LatticeLaw> latticeLaw(new LatticeLaw);
	latticeLaw->ensure2D   = false;
	latticeLaw->roughEdges = !smooth_edges_displacement;
	latticeLaw->calcTorsion= true;
	latticeLaw->tension_compression_different_stiffness = false;

	rootBody->engines.clear();
	rootBody->engines.push_back(boundDispatcher);
	rootBody->engines.push_back(latticeLaw);
	rootBody->engines.push_back(geometricalModelDispatcher);

	rootBody->initializers.clear();
	rootBody->initializers.push_back(boundDispatcher);
	rootBody->initializers.push_back(geometricalModelDispatcher);
}

void LatticeExampleCTData::positionRootBody(shared_ptr<Scene>& rootBody)
{
	rootBody->isDynamic             = false;

	Quaternionr q;
	q.FromAxisAngle( Vector3r(0,0,1),0);
	shared_ptr<LatticeSetParameters> physics(new LatticeSetParameters);
	physics->se3                    = Se3r(Vector3r(0,0,0),q);
	physics->beamGroupMask          = beamGroupMask;
	physics->nodeGroupMask          = nodeGroupMask;
	physics->useBendTensileSoftening= false;
	physics->useStiffnessSoftening  = false;

	shared_ptr<LatticeInteractingGeometry> set(new LatticeInteractingGeometry());


	set->diffuseColor		= Vector3r(0,0,1);

	shared_ptr<Aabb> aabb(new Aabb);
	aabb->diffuseColor		= Vector3r(0,0,1);

	shared_ptr<GeometricalModel> gm = YADE_PTR_CAST<GeometricalModel>(ClassFactory::instance().createShared("LatticeSetGeometry"));
	gm->diffuseColor 		= Vector3r(1,1,1);
	gm->wire 			= false;
	gm->shadowCaster 		= true;

	rootBody->shape	= YADE_PTR_CAST<Shape>(set);	
	rootBody->bound	= YADE_PTR_CAST<Bound>(aabb);
	rootBody->geometricalModel 	= gm;
	rootBody->physicalParameters 	= physics;
}

void LatticeExampleCTData::imposeTranslation(shared_ptr<Scene>& rootBody, Vector3r min, Vector3r max, Vector3r direction, Real displacement)
{
	shared_ptr<DisplacementEngine> translationCondition = shared_ptr<DisplacementEngine>(new DisplacementEngine);
	translationCondition->displacement  = displacement;
	direction.Normalize();
	translationCondition->translationAxis = direction;

	rootBody->engines.push_back(translationCondition);
	translationCondition->subscribedBodies.clear();

	BodyContainer::iterator bi    = rootBody->bodies->begin();
	BodyContainer::iterator biEnd = rootBody->bodies->end();
	for(  ; bi!=biEnd ; ++bi )
	{
		shared_ptr<Body> b = *bi;

		if( b->getGroupMask() & nodeGroupMask )
		{
			Vector3r pos = b->physicalParameters->se3.position;
			if(		   pos[0] > min[0] 
					&& pos[1] > min[1] 
					&& pos[2] > min[2] 
					&& pos[0] < max[0] 
					&& pos[1] < max[1] 
					&& pos[2] < max[2] 
					&& (b->getGroupMask() & nodeGroupMask)
			  )
			{
				b->isDynamic = false;
				b->geometricalModel->diffuseColor = Vector3r(2.0,2.0,0.0);
				translationCondition->subscribedBodies.push_back(b->getId());
			}
		}
	}
}

YADE_PLUGIN((LatticeExampleCTData));

YADE_REQUIRE_FEATURE(PHYSPAR);

