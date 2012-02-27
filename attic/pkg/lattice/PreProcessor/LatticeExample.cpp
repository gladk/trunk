/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "LatticeExample.hpp"

#include<yade/pkg-lattice/LatticeSetParameters.hpp>
#include<yade/pkg-lattice/LatticeBeamParameters.hpp>
#include<yade/pkg-lattice/LatticeBeamAngularSpring.hpp>
#include<yade/pkg-lattice/NonLocalDependency.hpp>
#include<yade/pkg-lattice/LatticeNodeParameters.hpp>
#include<yade/pkg-lattice/LineSegment.hpp>
#include<yade/pkg-lattice/LatticeLaw.hpp>
#include<yade/pkg-lattice/StrainRecorder.hpp>
#include<yade/pkg-lattice/NodeRecorder.hpp>
#include<yade/pkg-lattice/BeamRecorder.hpp>
#include<yade/pkg-lattice/MovingSupport.hpp>
#include<yade/pkg-lattice/MeasurePoisson.hpp>
#include<yade/pkg-lattice/NonLocalInitializer.hpp>
#include<yade/lib-base/yadeWm3Extra.hpp>

#include<yade/pkg-lattice/LatticeInteractingGeometry.hpp>

#include<yade/pkg-common/SphereModel.hpp>

#include<yade/core/Body.hpp>
#include<yade/core/BodyVector.hpp>
#include<yade/pkg-common/BoundDispatcher.hpp>
#include<yade/pkg-common/GeometricalModelMetaEngine.hpp>

#include<yade/pkg-common/Aabb.hpp>

#include<yade/core/BodyVector.hpp>
#include<yade/pkg-common/DisplacementEngine.hpp>
#include<yade/pkg-common/StateMetaEngine.hpp>
#include<yade/pkg-common/PhysicalActionApplier.hpp>

#include<yade/pkg-common/Quadrilateral.hpp>
#include<yade/pkg-common/ParticleParameters.hpp>

// Delaunay
#ifndef MINIWM3
	#include <Wm3Delaunay3.h>
	#include <Wm3Delaunay2.h>
	#include <Wm3Query.h>
#endif

YADE_REQUIRE_FEATURE(geometricalmodel)

using namespace boost;
using namespace std;


LatticeExample::LatticeExample() : FileGenerator()
{
        nodeGroupMask           = 1;
        beamGroupMask           = 2;
	quadGroupMask		= 4;
	fibreGroupMask		= 8;
	aggGroupMask		= 16;
	fibreBondGroupMask	= 32;
	bondGroupMask		= 64;
        
        speciemen_size_in_meters = Vector3r(0.1,0.1,0.0001);
        cellsizeUnit_in_meters   = 0.01;//0.003;
#ifndef MINIWM3
	use_Delaunay		 = true;
#else
	use_Delaunay		 = false;
#endif
        minAngle_betweenBeams_deg= 20.0;
        disorder_in_cellsizeUnit = Vector3r(0.6,0.6,0.0);
        maxLength_in_cellsizeUnit= 1.9;
        triangularBaseGrid       = true;
        triangularBaseGrid3D     = true;
        useNonLocalModel         = false;
        useBendTensileSoftening  = false;
        useStiffnessSoftening    = false;
        nonLocalL_in_cellsizeUnit= 1.1;         // l
                                
        crit_TensileStrain       = 1.0; // E_min
        crit_ComprStrain         = 0.5; // E_max
                        
        longitudalStiffness_noUnit= 1.0;        // k_l
        bendingStiffness_noUnit   = 0.6;        // k_b
        torsionalStiffness_noUnit = 0.6;        // k_t
        
        ignore_DOFs__better_is_OFF = true;
	ensure2D 		 = false;
        roughEdges 		 = false;
	calculate_Torsion_3D	 = false;
	quads			 = false;
        
	region_single_node_ABCDEF = false;

        region_A_min             = Vector3r(-0.006, 0.096,-1);
        region_A_max             = Vector3r( 0.16 , 0.16 , 1);
	direction_A 		 = Vector3r(0,1,0);
	blocked_xyz_A		 = Vector3r(0,1,1);
	displacement_A_meters	 = 0.0001;
	
	region_B_min 		 = Vector3r(-0.006,-0.006,-1);
	region_B_max 		 = Vector3r( 0.16 , 0.004, 1);
	direction_B 		 = Vector3r(0,-1,0);
	blocked_xyz_B		 = Vector3r(0,1,1);
	displacement_B_meters	 = 0.0001;

	region_C_min 		 = Vector3r(-0.006, 0.096,-1);
	region_C_max 		 = Vector3r( 0.16 , 0.16 ,-1);
	direction_C 		 = Vector3r(0,1,0);
	blocked_xyz_C		 = Vector3r(0,1,1);
	displacement_C_meters	 = 0.0001;
	
	region_D_min 		 = Vector3r(-0.006,-0.006,-1);
	region_D_max 		 = Vector3r( 0.16 , 0.004,-1);
        direction_D              = Vector3r(0,-1,0);
	blocked_xyz_D		 = Vector3r(0,1,1);
        displacement_D_meters    = 0.0001;
        
	region_E_min 		 = Vector3r(0,0,-1);
	region_E_max 		 = Vector3r(0,0,-1);
        direction_E              = Vector3r(0,1,0);
	blocked_xyz_E		 = Vector3r(0,1,1);
        displacement_E_meters    = 0.0;
        
	region_F_min 		 = Vector3r(0,0,-1);
	region_F_max 		 = Vector3r(0,0,-1);
        direction_F              = Vector3r(0,1,0);
	blocked_xyz_F		 = Vector3r(0,1,1);
        displacement_F_meters    = 0.0;
        
        strainRecorder_xz_plane  = -1;
        strainRecorder_node1     = Vector3r(0,0,0);
        strainRecorder_node2     = Vector3r(0,1,0);
        measurePoisson_node3     = Vector3r(0  ,0.1,0);
        measurePoisson_node4     = Vector3r(0.2,0.1,0);
        outputFile               = "./strains";
        poissonFile              = "./poisson";
        nodeRecorderFile         = "./npos";
        beamRecorderFile         = "./bstr";
        subscribedBodies.clear();
                
        regionDelete_A_min       = Vector3r(0,0,0);
        regionDelete_A_max       = Vector3r(0,0,0);
        regionDelete_B_min       = Vector3r(0,0,0);
        regionDelete_B_max       = Vector3r(0,0,0);
        regionDelete_C_min       = Vector3r(0,0,0);
        regionDelete_C_max       = Vector3r(0,0,0);
        regionDelete_D_min       = Vector3r(0,0,0);
        regionDelete_D_max       = Vector3r(0,0,0);
        regionDelete_E_min       = Vector3r(0,0,0);
        regionDelete_E_max       = Vector3r(0,0,0);
        regionDelete_F_min       = Vector3r(0,0,0);
        regionDelete_F_max       = Vector3r(0,0,0);

        regionDelete_1_min       = Vector3r(0,0,0);
        regionDelete_1_max       = Vector3r(0,0,0);
        regionDelete_2_min       = Vector3r(0,0,0);
        regionDelete_2_max       = Vector3r(0,0,0);
        regionDelete_3_min       = Vector3r(0,0,0);
        regionDelete_3_max       = Vector3r(0,0,0);
        regionDelete_4_min       = Vector3r(0,0,0);
        regionDelete_4_max       = Vector3r(0,0,0);
        regionDelete_5_min       = Vector3r(0,0,0);
        regionDelete_5_max       = Vector3r(0,0,0);

        nonDestroy_A_min         = Vector3r(0,0,0);
        nonDestroy_A_max         = Vector3r(0,0,0);
        nonDestroy_B_min         = Vector3r(0,0,0);
        nonDestroy_B_max         = Vector3r(0,0,0);
	nonDestroy_stiffness	 = 10.0;

	CT_A_min		 = Vector3r(0,0,-1);
	CT_A_max		 = Vector3r(0,0,-1);
	CT_B_min		 = Vector3r(0,0,-1);
	CT_B_max		 = Vector3r(0,0,-1);
	CT			 = 1;

        useAggregates            = false;
	no_Agg_outside		 = false;
        aggregatePercent         = 40;
        aggregateMeanDiameter    = cellsizeUnit_in_meters*1;
        aggregateSigmaDiameter   = cellsizeUnit_in_meters*2;
        aggregateMinDiameter     = cellsizeUnit_in_meters*2;
        aggregateMaxDiameter     = cellsizeUnit_in_meters*4;
        // MaterialParameters of aggregate
        agg_longStiffness_noUnit = 3.0;                                 // k_l aggregate
        agg_bendStiffness_noUnit = 2.1;                                 // k_b aggregate
        agg_torsStiffness_noUnit = 2.1;                                 // k_b aggregate
        agg_critCompressStrain   = 100.0;                               // E.c aggregate
        agg_critTensileStrain    = 50.0;                                // E.l aggregate
        // MaterialParameters of bond
        bond_longStiffness_noUnit= 0.7;                                 // k_l bond
        bond_bendStiffness_noUnit= 0.28;                                // k_b bond
        bond_torsStiffness_noUnit= 0.28;                                // k_b bond
        bond_critCompressStrain  = 100.0;                               // E.c bond
        bond_critTensileStrain   = 50.0;                                // E.l bond
	// MaterialParameters of steel fibres
        fibre_longStiffness_noUnit=8.0;      // k_l fibre
        fibre_bendStiffness_noUnit=5.6;      // k_b fibre
        fibre_torsStiffness_noUnit=5.6;      // k_t fibre
        fibre_critCompressStrain = 1;        // E.c fibre
        fibre_critTensileStrain  = 0.0055;   // E.l fibre
	fibre_count		 = 0;
	beams_per_fibre		 = 10;
	fibre_allows		 = 0.5;
	//fibre_irregularity_noUnit= 5;
	fibre_balancing_iterations= 300;
        // MaterialParameters of fibre bond
        fibre_bond_longStiffness_noUnit= 0.7;         // k_l fibre bond
        fibre_bond_bendStiffness_noUnit= 0.28;        // k_b fibre bond
        fibre_bond_torsStiffness_noUnit= 0.28;        // k_t fibre bond
        fibre_bond_critCompressStrain  = 100.0;       // E.c fibre bond
        fibre_bond_critTensileStrain   = 50.0;        // E.l fibre bond
	
	record_only_matrix = false;
	nodeRec_A_min=Vector3r(0,0,0);
	nodeRec_A_max=Vector3r(0,0,0);
	nodeRec_B_min=Vector3r(0,0,0);
	nodeRec_B_max=Vector3r(0,0,0);
	nodeRec_C_min=Vector3r(0,0,0);
	nodeRec_C_max=Vector3r(0,0,0);
	nodeRec_D_min=Vector3r(0,0,0);
	nodeRec_D_max=Vector3r(0,0,0);
	nodeRec_E_min=Vector3r(0,0,0);
	nodeRec_E_max=Vector3r(0,0,0);
	nodeRec_F_min=Vector3r(0,0,0);
	nodeRec_F_max=Vector3r(0,0,0);
	
	beamRec_A_pos=Vector3r(0,0,0);
	beamRec_B_pos=Vector3r(0,0,0);
	beamRec_C_pos=Vector3r(0,0,0);
	beamRec_D_pos=Vector3r(0,0,0);
	beamRec_E_pos=Vector3r(0,0,0);

	beamRec_A_range=0;
	beamRec_B_range=0;
	beamRec_C_range=0;
	beamRec_D_range=0;
	beamRec_E_range=0;

	beamRec_A_dir=-1;
	beamRec_B_dir=-1;
	beamRec_C_dir=-1;
	beamRec_D_dir=-1;
	beamRec_E_dir=-1;
		
	movSupp_A_pos=Vector3r(0,0,0);
	movSupp_A_range=0;
	movSupp_A_dir=-1;

	movSupp_B_pos=Vector3r(0,0,0);
	movSupp_B_range=0;
	movSupp_B_dir=-1;

	movSupp_C_pos=Vector3r(0,0,0);
	movSupp_C_range=0;
	movSupp_C_dir=-1;

	movSupp_D_pos=Vector3r(0,0,0);
	movSupp_D_range=0;
	movSupp_D_dir=-1;

	fibres_horizontal = false;
	fibres_vertical = false;
}


LatticeExample::~LatticeExample()
{

}



bool LatticeExample::generate()
{
	no_Agg_outside = false;
	std::cerr << "no_Agg_outside parameter is not used.\n";

	fibres_total=0;
	matrix_total=0;
	beam_total=0;
	rootBody = shared_ptr<Scene>(new Scene);
	createActors(rootBody);
	positionRootBody(rootBody);
	
	
	shared_ptr<Body> body;
	
	Vector3r nbNodes = speciemen_size_in_meters / cellsizeUnit_in_meters;
	if(triangularBaseGrid)
		nbNodes[1] *= 1.15471; // bigger by sqrt(3)/2 factor
	if(triangularBaseGrid3D)
		nbNodes[2] *= 1.22475; // bigger by (1/3)*(sqrt(6)) factor


	bool FLAT = speciemen_size_in_meters[2]<cellsizeUnit_in_meters;
	std::vector<Vector2r> vert2; vert2.clear(); // Delaunay
	std::vector<Vector3r> vert3; vert3.clear(); // Delaunay

        unsigned int totalNodesCount = 0;
	{
		if(fibre_count > 0.0)
			makeFibres();

		setStatus("creating nodes...");
		float all = nbNodes[0]*nbNodes[1]*nbNodes[2];
		float current = 0.0;

		for( int j=0 ; j<=nbNodes[1] ; j++ )
		{
			if(shouldTerminate()) return false;

			for( int i=0 ; i<=nbNodes[0] ; i++ )
				for( int k=0 ; k<=nbNodes[2] ; k++)
				{
					setProgress(current++/all);

					shared_ptr<Body> node;
					if(createNode(node,i,j,k) || quads)
					{
						if( fibre_count < 0.5 || fibreAllows(node->physicalParameters->se3.position))
						{
						rootBody->bodies->insert(node), ++totalNodesCount;
						if(FLAT)
							vert2.push_back(Vector2r(node->physicalParameters->se3.position[0],node->physicalParameters->se3.position[1])); // Delaunay
						else
							vert3.push_back(node->physicalParameters->se3.position); // Delaunay
						}
					}
				}
		}
		if(fibre_count > 0.0)
		{
			setStatus("creating fibre nodes...");
			for(int i = 0 ; i < fibre_count ; ++i)
		        {
				Vector3r pos = fibres[i].first;
				Vector3r del = fibres[i].second;
				for(int j = 0 ; j < beams_per_fibre ; ++j)
				{
					shared_ptr<Body> node;
					if(createNodeXYZ(
						node,
						pos[0] + 1.0*j*del[0],
						pos[1] + 1.0*j*del[1],
						pos[2] + 1.0*j*del[2] ))
					{
						rootBody->bodies->insert(node), ++totalNodesCount;
						if(FLAT)
							vert2.push_back(Vector2r(node->physicalParameters->se3.position[0],node->physicalParameters->se3.position[1])); // Delaunay
						else
							vert3.push_back(node->physicalParameters->se3.position); // Delaunay
					}
				}
			}
		}

	}

        BodyVector bc;
        bc.clear();

	BodyContainer::iterator bi    = rootBody->bodies->begin();
	BodyContainer::iterator bi2;
	BodyContainer::iterator biEnd = rootBody->bodies->end();

#ifndef MINIWM3
	if(use_Delaunay) // create beams, Delaunay
	{
		std::set< std::pair<int,int> > pairs;
		int I,J;
		if(FLAT)
		{
			setProgress(0); setStatus("Delaunay 2d...");
			Delaunay2<Real> del2(vert2.size(),&(vert2[0]),cellsizeUnit_in_meters/50,false, Query::QT_INTEGER); // Delaunay

			int del2_i = 0;
			int del2_ind[3];
			while(del2.GetIndexSet(del2_i++,del2_ind))
			{
				if(shouldTerminate()) return false;
				setProgress((float)del2_i/(float)del2.GetSimplexQuantity());
				shared_ptr<Body> beam;

				for(int zzz=0 ; zzz<3 ; ++zzz)
				{
					I=del2_ind[zzz];J=del2_ind[(zzz+1)%3];
					if(J>I) std::swap(I,J);
					if(pairs.insert(std::make_pair(I,J)).second)
					{
						createBeam(beam,I,J);
						if(calcBeamPositionOrientationLength(beam)<maxLength_in_cellsizeUnit*cellsizeUnit_in_meters)
							if(notDeleted(beam->physicalParameters->se3.position)) 
								bc.insert(beam);
					}
				}
			}
			setStatus("Delaunay 2d... done");
		}
		else
		{
			setProgress(0); setStatus("Delaunay 3d...");
			Delaunay3<Real> del3(vert3.size(),&(vert3[0]),cellsizeUnit_in_meters/50,false, Query::QT_INTEGER); // Delaunay
			
			int del3_i = 0;
			int del3_ind[4];
			while(del3.GetIndexSet(del3_i++,del3_ind))
			{
				if(shouldTerminate()) return false;
				setProgress((float)del3_i/(float)del3.GetSimplexQuantity());
				shared_ptr<Body> beam;

				for(int zzz=0 ; zzz<4 ; ++zzz)
				{
					I=del3_ind[zzz];J=del3_ind[(zzz+1)%4];
					if(J>I) std::swap(I,J);
					if(pairs.insert(std::make_pair(I,J)).second)
					{
						createBeam(beam,I,J);
						if(calcBeamPositionOrientationLength(beam)<maxLength_in_cellsizeUnit*cellsizeUnit_in_meters)
							if(notDeleted(beam->physicalParameters->se3.position)) 
								bc.insert(beam);
					}
				}
			}
			setStatus("Delaunay 3d... done");
		}

	}
#else
	if(use_Delaunay) // create beams, Delaunay not avauilable
	{
		message = "ERROR:\n\nYade was compiled without full installation of wildmagic-dev (Wm3 foundation library), can't use Delaunay. Please disable option use_Delaunay or compile yade with full wm3 library installed. Don't forget CPPPATH=/usr/include/wm3 scons parameter.";
		return false;
	}
#endif
	else
	{ //  create beams, old method
	int beam_counter = 0;
	float nodes_a=0;
	float nodes_all = rootBody->bodies->size();
	setStatus("creating beams...");
	double r  = maxLength_in_cellsizeUnit*cellsizeUnit_in_meters;
	double r2 = std::pow(r,2);
	for(  ; bi!=biEnd ; ++bi )  // loop over all nodes, to create beams
	{
		Body* bodyA = (*bi).get(); // first_node
	
		bi2 = bi;
		++bi2;
		nodes_a+=1.0;
		
		for( ; bi2!=biEnd ; ++bi2 )
		{
			if(shouldTerminate()) return false;

			Body* bodyB = (*bi2).get(); // all other nodes
			// warning - I'm assuming that there are ONLY Nodes in the rootBody
			LatticeNodeParameters* a = static_cast<LatticeNodeParameters*>(bodyA->physicalParameters.get());
			LatticeNodeParameters* b = static_cast<LatticeNodeParameters*>(bodyB->physicalParameters.get());
			
			//if ((a->se3.position - b->se3.position).SquaredLength() < std::pow(maxLength_in_cellsizeUnit*cellsizeUnit_in_meters,2) )  

			if (	   ( std::abs(a->se3.position[0] - b->se3.position[0])<= r )
				&& ( std::abs(a->se3.position[1] - b->se3.position[1])<= r )
				&& ( std::abs(a->se3.position[2] - b->se3.position[2])<= r )
				&& ((a->se3.position - b->se3.position).SquaredLength() < r2 ) )  
			{
				shared_ptr<Body> beam;
				createBeam(beam,bodyA->getId(),bodyB->getId());
				calcBeamPositionOrientationLength(beam);
				if(checkMinimumAngle(bc,beam))
				{
					if( ++beam_counter % 100 == 0 )
					{
						cerr << "creating beam: " << beam_counter << " , " << ((nodes_a/nodes_all)*100.0)  << " %\n"; 
						setProgress(std::pow(nodes_a/nodes_all,2));
					}
					
					if(notDeleted(beam->physicalParameters->se3.position)) 
						bc.insert(beam);
				}
			}
                }
        }
	} // beams are created


	setStatus("strain recorder nodes...");
        { // subscribe two nodes, that are monitored by strain recoder to get a measure of length
                bi    = rootBody->bodies->begin();
                biEnd = rootBody->bodies->end();
                int node1Id=-1; Real len1=100000; Real tmpLen;
                int node2Id=-1; Real len2=100000;
                for(  ; bi!=biEnd ; ++bi )  // loop over all nodes, to find those closest to strainRecorder_node1 and strainRecorder_node2
                {
			if(shouldTerminate()) return false;

                        Body* body = (*bi).get();
                        LatticeNodeParameters* node = YADE_CAST<LatticeNodeParameters*>(body->physicalParameters.get());
                        tmpLen = ( strainRecorder_node1 - node->se3.position ).SquaredLength();
                        if(tmpLen < len1) len1=tmpLen, node1Id=body->getId();
                        tmpLen = ( strainRecorder_node2 - node->se3.position ).SquaredLength();
                        if(tmpLen < len2) len2=tmpLen, node2Id=body->getId();
                }
                assert(node1Id!=-1 && node2Id!=-1);
                subscribedBodies.push_back(node1Id); // bottom
                subscribedBodies.push_back(node2Id); // upper
                
                measurePoisson->bottom =node1Id; // bottom
                measurePoisson->upper  =node2Id; // upper
        }
        
	setStatus("poisson recorder nodes...");
        { // subscribe two nodes, that are monitored by MeasurePoisson to get a measure of poisson
                bi    = rootBody->bodies->begin();
                biEnd = rootBody->bodies->end();
                int node3Id=-1; Real len1=100000; Real tmpLen;
                int node4Id=-1; Real len2=100000;
                for(  ; bi!=biEnd ; ++bi )  // loop over all nodes, to find those closest to measurePoisson_node3 and measurePoisson_node4
                {
			if(shouldTerminate()) return false;

                        Body* body = (*bi).get();
                        LatticeNodeParameters* node = YADE_CAST<LatticeNodeParameters*>(body->physicalParameters.get());
                        tmpLen = ( measurePoisson_node3 - node->se3.position ).SquaredLength();
                        if(tmpLen < len1) len1=tmpLen, node3Id=body->getId();
                        tmpLen = ( measurePoisson_node4 - node->se3.position ).SquaredLength();
                        if(tmpLen < len2) len2=tmpLen, node4Id=body->getId();
                }
                assert(node3Id!=-1 && node4Id!=-1);
                measurePoisson->left   =node3Id; // left
                measurePoisson->right  =node4Id; // right
                measurePoisson->horizontal              =   (*(rootBody->bodies))[node4Id]->physicalParameters->se3.position[0] 
                                                          - (*(rootBody->bodies))[node3Id]->physicalParameters->se3.position[0];
        }


        bi    = bc.begin();
        biEnd = bc.end();
	float all_bodies = bc.size();
	setStatus("Copy beam data...");
        for( int curr=0 ; bi!=biEnd ; ++bi,++curr )  // loop over all newly created beams ...
        {
		if(curr % 100 == 0)
		{
			if(shouldTerminate()) return false;
			setProgress((float)(curr)/(float)(all_bodies));
		}

                shared_ptr<Body> b = *bi;
                rootBody->bodies->insert(b); // .. to insert them into rootBody
                
                // attach them to strain recorder
                
                Real xz_plane = strainRecorder_xz_plane;
                Real pos1 = (*(rootBody->bodies))[static_cast<LatticeBeamParameters*>( b->physicalParameters.get() )->id1]->physicalParameters->se3.position[1]; // beam first node
                Real pos2 = (*(rootBody->bodies))[static_cast<LatticeBeamParameters*>( b->physicalParameters.get() )->id2]->physicalParameters->se3.position[1]; // beam second node
                if( pos1 < pos2 )
                        std::swap(pos1,pos2); // make sure that pos1 is bigger 
                if(        pos1 > xz_plane
                        && pos2 < xz_plane 
                        )
                        subscribedBodies.push_back( b->getId() ); // beam crosses the plane!
        }
        strainRecorder->subscribedBodies        = subscribedBodies;

        strainRecorder->initialLength           =   (*(rootBody->bodies))[subscribedBodies[0]]->physicalParameters->se3.position[1] 
                                                  - (*(rootBody->bodies))[subscribedBodies[1]]->physicalParameters->se3.position[1];
        
        measurePoisson->vertical                =   (*(rootBody->bodies))[subscribedBodies[1]]->physicalParameters->se3.position[1] 
                                                  - (*(rootBody->bodies))[subscribedBodies[0]]->physicalParameters->se3.position[1];
                
	setStatus("Find neighbours...");
        { // remember what node is in contact with what beams
                //    node                   beams
		connections.resize(totalNodesCount);
		
		bi    = rootBody->bodies->begin();
		biEnd = rootBody->bodies->end();
		float all_bodies = rootBody->bodies->size();
		
		for( int curr=0 ; bi!=biEnd ; ++bi,++curr )  // loop over all beams
		{
			if(curr % 100 == 0)
			{
				if(shouldTerminate()) return false;
				setProgress((float)(curr)/(float)(all_bodies));
			}

			Body* body = (*bi).get();
			if( ! ( body->getGroupMask() & beamGroupMask ) )
				continue; // skip non-beams
			
			LatticeBeamParameters* beam = static_cast<LatticeBeamParameters*>(body->physicalParameters.get() );
			connections[beam->id1].push_back(body->getId());
			connections[beam->id2].push_back(body->getId());
		}
	}

	setStatus("angular springs...");
	{ // create angular springs between beams
		bi    = rootBody->bodies->begin();
		biEnd = rootBody->bodies->end();
		float all_bodies = rootBody->bodies->size();
		for( int curr=0 ; bi!=biEnd ; ++bi,++curr )  // loop over all beams
		{
			if(curr % 100 == 0)
			{
				if(shouldTerminate()) return false;
				setProgress((float)(curr)/(float)(all_bodies));
			}

			Body* body = (*bi).get();
			if( ! ( body->getGroupMask() & beamGroupMask ) )
				continue; // skip non-beams
				
			calcBeamAngles(body,rootBody->bodies.get(),rootBody->interactions.get());
                }
        }
	
	if(quads)
	{
		float all = (nbNodes[0]-1)*(nbNodes[1]-1)*(nbNodes[2]-1);
		float current = 0.0;
		setStatus("quadrilaterals...");
		for( int j=1 ; j<=nbNodes[1] ; j++ )
		{
			if(shouldTerminate()) return false;

			for( int i=1 ; i<=nbNodes[0] ; i++ )
			{
				shared_ptr<Body> quad;
				if(createQuad(quad,i,j,nbNodes))
					rootBody->bodies->insert(quad);
				setProgress(current++/all);
			}
		}
	};
        
        imposeTranslation(rootBody,region_A_min,region_A_max,direction_A,displacement_A_meters,blocked_xyz_A);
        imposeTranslation(rootBody,region_B_min,region_B_max,direction_B,displacement_B_meters,blocked_xyz_B);
        imposeTranslation(rootBody,region_C_min,region_C_max,direction_C,displacement_C_meters,blocked_xyz_C);
        imposeTranslation(rootBody,region_D_min,region_D_max,direction_D,displacement_D_meters,blocked_xyz_D);
        imposeTranslation(rootBody,region_E_min,region_E_max,direction_E,displacement_E_meters,blocked_xyz_E);
        imposeTranslation(rootBody,region_F_min,region_F_max,direction_F,displacement_F_meters,blocked_xyz_F);

	beam_total=bc.size();
        if(useAggregates) addAggregates(rootBody);
	setStatus("making fibres...");
	if(fibre_count > 0.0) makeFibreBeams(rootBody);
	if(shouldTerminate()) return false;
        
        nonDestroy(rootBody,nonDestroy_A_min,nonDestroy_A_max);
        nonDestroy(rootBody,nonDestroy_B_min,nonDestroy_B_max);
        
	modifyCT(rootBody,CT_A_min,CT_A_max);
        modifyCT(rootBody,CT_B_min,CT_B_max);

        message="Number of nodes created:\n" + lexical_cast<string>(nbNodes[0]) + ","
                                            + lexical_cast<string>(nbNodes[1]) + ","
                                            + lexical_cast<string>(nbNodes[2]) + ",\n"
	     + "Number of beams: " + lexical_cast<string>(bc.size()) + "\n"
	     + "Fibres total: " + lexical_cast<string>(fibres_total) + "\n"
	     + "Matrix beams total: " + lexical_cast<string>(matrix_total) + "\n"
//	     + "Fibres/matrix %: " + lexical_cast<string>(100.0*fibres_total/matrix_total) + "\n"
	     + "Fibres/all other beams %: " + lexical_cast<string>(100.0*fibres_total/(beam_total-fibres_total)) + "\n"
	     + "\nNOTE: sometimes it can look better when 'drawWireFrame' is enabled in Display tab.";

        cerr << "finished.. saving\n" << message << "\n";
	
	return true;
}

/// returns true if angle is bigger than minAngle_betweenBeams_deg
bool LatticeExample::checkAngle( Vector3r a, Vector3r& b)
{
	Quaternionr al;
	al.Align(a,b);
	Vector3r axis;
	Real angle;
	al.ToAxisAngle(axis, angle);
	angle *= 180.0/Mathr::PI ;
//	cerr << " angle: " << angle << "\n";
	return angle > minAngle_betweenBeams_deg;
}

/// returns true if angle is bigger than minAngle_betweenBeams_deg
bool LatticeExample::checkMinimumAngle(BodyVector& bc,shared_ptr<Body>& body)
{
	bool answer = true;
	LatticeBeamParameters* newBeam = static_cast<LatticeBeamParameters*>(body->physicalParameters.get());
	
	BodyContainer::iterator bi    = bc.begin();
	BodyContainer::iterator biEnd = bc.end();
	for(  ; bi!=biEnd ; ++bi )  // loop over all beams - they MUST be beams, for static_cast<> 
	{ 
		LatticeBeamParameters* oldBeam = static_cast<LatticeBeamParameters*>((*bi)->physicalParameters.get());
		if( 	   (oldBeam->id1 == newBeam->id1)
			|| (oldBeam->id2 == newBeam->id2))
			answer = answer && checkAngle(   oldBeam->direction ,  newBeam->direction );
		if( 	   (oldBeam->id2 == newBeam->id1)
			|| (oldBeam->id1 == newBeam->id2))
			answer = answer && checkAngle( - oldBeam->direction ,  newBeam->direction );
	} 
	return answer;
}

bool LatticeExample::createNodeXYZ(shared_ptr<Body>& body, Real x, Real y, Real z)
{
	body = shared_ptr<Body>(new Body(body_id_t(0),nodeGroupMask));
	shared_ptr<LatticeNodeParameters> physics(new LatticeNodeParameters);
	shared_ptr<SphereModel> gSphere(new SphereModel);
	
	Quaternionr q;
	q.FromAxisAngle( Vector3r(Mathr::UnitRandom(),Mathr::UnitRandom(),Mathr::UnitRandom()) , Mathr::UnitRandom()*Mathr::PI );
	
	Vector3r position(x,y,z);

	Real radius 			= cellsizeUnit_in_meters*0.05;
	
	body->isDynamic			= true;
	
	physics->se3			= Se3r(position,q);

	gSphere->radius			= radius;
	gSphere->diffuseColor		= Vector3r(1.8,1.8,0.0);
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

bool LatticeExample::createNode(shared_ptr<Body>& body, int i, int j, int k)
{
	body = shared_ptr<Body>(new Body(body_id_t(0),nodeGroupMask));
	shared_ptr<LatticeNodeParameters> physics(new LatticeNodeParameters);
	shared_ptr<SphereModel> gSphere(new SphereModel);
	
	Quaternionr q;
	q.FromAxisAngle( Vector3r(Mathr::UnitRandom(),Mathr::UnitRandom(),Mathr::UnitRandom()) , Mathr::UnitRandom()*Mathr::PI );
	
	float  triang_x = triangularBaseGrid   ? ((static_cast<float>(j%2))*0.5)+i : i+0;
	double triang_y = triangularBaseGrid   ? 0.86602540378443864676*j          : j*1; // sqrt(3)/2
	double triang_z = k;

	if(triangularBaseGrid3D)
	{
		triang_z *= 0.81649658092772603273;
		switch(k%3)
		{
			case 0 : triang_x += 0.5; triang_y +=     0.86602540378443864676/3.0; break;
			case 1 :                  triang_y += 2.0*0.86602540378443864676/3.0; break;
			case 2 : break;
		}
	}
	
	Vector3r position		= ( Vector3r(triang_x,triang_y,triang_z)
					  + Vector3r( 	  Mathr::SymmetricRandom()*disorder_in_cellsizeUnit[0]
					  		, Mathr::SymmetricRandom()*disorder_in_cellsizeUnit[1]
							, Mathr::SymmetricRandom()*disorder_in_cellsizeUnit[2]
						    ) * 0.5 // *0.5 because symmetricRandom is (-1,1), and disorder is whole size where nodes can appear
					  )*cellsizeUnit_in_meters;

	Real radius 			= cellsizeUnit_in_meters*0.05;
	
	body->isDynamic			= true;
	
	physics->se3			= Se3r(position,q);

	gSphere->radius			= radius;
	gSphere->diffuseColor		= Vector3r(0.8,0.8,0.8);
	gSphere->wire			= false;
	gSphere->shadowCaster		= false;
	
        body->geometricalModel          = gSphere;
        body->physicalParameters        = physics;

//////////////////////
// some custom coordinates for testing: - FIXME - it will be useful to be able to get coordinates from outside
//static int zzzzz=0;
//switch(zzzzz%5)
//{
//      case 0  : physics->se3.position=Vector3r(0,0,0); break;
//      case 1  : physics->se3.position=Vector3r(3,0,0); break;
//      case 2  : physics->se3.position=Vector3r(8,0,0); break;
//      case 3  : physics->se3.position=Vector3r(0,4,0); break;
//      case 4  : physics->se3.position=Vector3r(3,4,0); break;
//};
//zzzzz++;
//
////      case 0  : physics->se3.position=Vector3r(0.4,1.5,0); break;
////      case 1  : physics->se3.position=Vector3r(0.8,0.6,0); break;
////      case 2  : physics->se3.position=Vector3r(0  ,0  ,0); break;
////      case 3  : physics->se3.position=Vector3r(1.6,0.5,0); break;
////      case 4  : physics->se3.position=Vector3r(2.0,0  ,0); break;
////
////////////////////////
        
	if( 	   position[0] >= speciemen_size_in_meters[0] 
		|| position[1] >= speciemen_size_in_meters[1]
		|| position[2] >= speciemen_size_in_meters[2] )
		return false;

        return true;
}


bool LatticeExample::createQuad(shared_ptr<Body>& quad, int i, int j, Vector3r nbNodes)
{
	quad = shared_ptr<Body>(new Body(body_id_t(0),quadGroupMask));
	shared_ptr<Quadrilateral> gQuad(new Quadrilateral( (j-1)*(int)std::floor(nbNodes[0]+1)+i-1 , j*(int)std::floor(nbNodes[0]+1)+i-1 , j*(int)std::floor(nbNodes[0]+1)+i , (j-1)*(int)std::floor(nbNodes[0]+1)+i, rootBody.get()));
	// FIXME FIXME - this is an empty class. not needed at all.
	// all this Quadrilateral stuff in fact does not fit current design at all.
	shared_ptr<ParticleParameters> physics(new ParticleParameters);
	
	quad->isDynamic			= false;

	physics->se3			= Se3r(Vector3r(0,0,0),Quaternionr(1,0,0,0));
	physics->mass			= 0;
	physics->velocity		= Vector3r(0,0,0);

	gQuad->diffuseColor		= Vector3r(0.0,0.0,0.0);
	gQuad->wire			= false;
	gQuad->shadowCaster		= false;
	
        quad->geometricalModel          = gQuad;
	quad->physicalParameters	= physics;

        return true;
}


void LatticeExample::createBeam(shared_ptr<Body>& body, int i, int j)
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


Real LatticeExample::calcBeamPositionOrientationLength(shared_ptr<Body>& body)
{
//	LatticeBeamParameters* beam = static_cast<LatticeBeamParameters*>(body->physicalParameters.get());

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
        
        beam->criticalTensileStrain     = crit_TensileStrain;
        beam->criticalCompressiveStrain = crit_ComprStrain;
        beam->longitudalStiffness       = longitudalStiffness_noUnit;
        beam->bendingStiffness          = bendingStiffness_noUnit;
        beam->torsionalStiffness        = torsionalStiffness_noUnit;
        
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

void LatticeExample::calcAxisAngle(LatticeBeamParameters* beam1, BodyContainer* bodies, int otherId, InteractionContainer* ints, int thisId)
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

		// FIXME
		if( dir1.Dot(beam1->otherDirection) < 0.999999 )
			angularSpring->initialOffPlaneAngle1   *= -1.0;//, angularSpring->offPlaneSwap1 = true;
		if( dir2.Dot(beam2->otherDirection) < 0.999999 )
			angularSpring->initialOffPlaneAngle2   *= -1.0;//, angularSpring->offPlaneSwap2 = true;
	
		interaction->interactionPhysics 	= angularSpring;
		ints->insert(interaction);
	}
}

void LatticeExample::calcBeamAngles(Body* body, BodyContainer* bodies, InteractionContainer* ints)
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

void LatticeExample::createActors(shared_ptr<Scene>& )
{
        shared_ptr<BoundDispatcher> boundDispatcher   = shared_ptr<BoundDispatcher>(new BoundDispatcher);

        shared_ptr<GeometricalModelMetaEngine> geometricalModelDispatcher       = shared_ptr<GeometricalModelMetaEngine>(new GeometricalModelMetaEngine);
        geometricalModelDispatcher->add("LatticeSet2LatticeBeams");
        
        strainRecorder = shared_ptr<StrainRecorder>(new StrainRecorder);
        strainRecorder->outputFile              = outputFile;
        strainRecorder->interval                = 10;
        
        measurePoisson = shared_ptr<MeasurePoisson>(new MeasurePoisson);
        measurePoisson->outputFile              = poissonFile;
        measurePoisson->interval                = 10;
        
/*
NodeRecorder nnnn;
	std::string	 outputFile;
	unsigned int	 interval;
	std::list<std::pair<Vector3r,Vector3r> > regions; // a list of min/max pairs describing each region.
BeamRecorder bbbb;
	std::string	 outputFile;
	unsigned int	 interval;
	std::list<std::pair<Vector3r,std::pair< Real, int > > > sections; // a list of section planes: <a midpoint, half length, direction >
*/
        nodeRecorder   = shared_ptr<NodeRecorder>(new NodeRecorder);
        nodeRecorder->outputFile                = nodeRecorderFile;
        nodeRecorder->interval                  = 10;
	if(record_only_matrix)
		nodeRecorder->only_this_stiffness	= longitudalStiffness_noUnit;
	else
		nodeRecorder->only_this_stiffness	= -1;
	if(nodeRec_A_min != nodeRec_A_max) nodeRecorder->regions.push_back(std::make_pair(nodeRec_A_min,nodeRec_A_max));
	if(nodeRec_B_min != nodeRec_B_max) nodeRecorder->regions.push_back(std::make_pair(nodeRec_B_min,nodeRec_B_max));
	if(nodeRec_C_min != nodeRec_C_max) nodeRecorder->regions.push_back(std::make_pair(nodeRec_C_min,nodeRec_C_max));
	if(nodeRec_D_min != nodeRec_D_max) nodeRecorder->regions.push_back(std::make_pair(nodeRec_D_min,nodeRec_D_max));
	if(nodeRec_E_min != nodeRec_E_max) nodeRecorder->regions.push_back(std::make_pair(nodeRec_E_min,nodeRec_E_max));
	if(nodeRec_F_min != nodeRec_F_max) nodeRecorder->regions.push_back(std::make_pair(nodeRec_F_min,nodeRec_F_max));
        
        beamRecorder = shared_ptr<BeamRecorder>(new BeamRecorder);
        beamRecorder->outputFile                = beamRecorderFile;
        beamRecorder->interval                  = 10;
	if(beamRec_A_dir != -1) beamRecorder->sections.push_back(std::make_pair(beamRec_A_pos,std::make_pair(beamRec_A_range,beamRec_A_dir)));
	if(beamRec_B_dir != -1) beamRecorder->sections.push_back(std::make_pair(beamRec_B_pos,std::make_pair(beamRec_B_range,beamRec_B_dir)));
	if(beamRec_C_dir != -1) beamRecorder->sections.push_back(std::make_pair(beamRec_C_pos,std::make_pair(beamRec_C_range,beamRec_C_dir)));
	if(beamRec_D_dir != -1) beamRecorder->sections.push_back(std::make_pair(beamRec_D_pos,std::make_pair(beamRec_D_range,beamRec_D_dir)));
	if(beamRec_E_dir != -1) beamRecorder->sections.push_back(std::make_pair(beamRec_E_pos,std::make_pair(beamRec_E_range,beamRec_E_dir)));
        
        movingSupport   = shared_ptr<MovingSupport>(new MovingSupport);
	if(movSupp_A_dir != -1) movingSupport->sections.push_back(std::make_pair(movSupp_A_pos,std::make_pair(movSupp_A_dir,movSupp_A_range)));
	if(movSupp_B_dir != -1) movingSupport->sections.push_back(std::make_pair(movSupp_B_pos,std::make_pair(movSupp_B_dir,movSupp_B_range)));
	if(movSupp_C_dir != -1) movingSupport->sections.push_back(std::make_pair(movSupp_C_pos,std::make_pair(movSupp_C_dir,movSupp_C_range)));
	if(movSupp_D_dir != -1) movingSupport->sections.push_back(std::make_pair(movSupp_D_pos,std::make_pair(movSupp_D_dir,movSupp_D_range)));
 	
	shared_ptr<LatticeLaw> latticeLaw(new LatticeLaw);

        latticeLaw->ensure2D   = ensure2D;
        latticeLaw->roughEdges = roughEdges;
	latticeLaw->backward_compatible = ignore_DOFs__better_is_OFF;


        latticeLaw->calcTorsion= calculate_Torsion_3D;
	latticeLaw->tension_compression_different_stiffness = true;
	latticeLaw->respect_non_destroy = nonDestroy_stiffness;
        
        rootBody->engines.clear();
        rootBody->engines.push_back(boundDispatcher);
        rootBody->engines.push_back(latticeLaw);
        rootBody->engines.push_back(geometricalModelDispatcher);
        rootBody->engines.push_back(strainRecorder);
        rootBody->engines.push_back(measurePoisson);
        rootBody->engines.push_back(nodeRecorder);
        rootBody->engines.push_back(beamRecorder);
	rootBody->engines.push_back(movingSupport);
        
        rootBody->initializers.clear();
        rootBody->initializers.push_back(boundDispatcher);
        rootBody->initializers.push_back(geometricalModelDispatcher);

        if(useNonLocalModel)
        {
                shared_ptr<NonLocalInitializer> nonLocalInitializer(new NonLocalInitializer);
                nonLocalInitializer->range = nonLocalL_in_cellsizeUnit * cellsizeUnit_in_meters;
                rootBody->initializers.push_back(nonLocalInitializer);
        }
}       

void LatticeExample::positionRootBody(shared_ptr<Scene>& rootBody)
{
        rootBody->isDynamic             = false;

	Quaternionr q;
	q.FromAxisAngle( Vector3r(0,0,1),0);
	shared_ptr<LatticeSetParameters> physics(new LatticeSetParameters);
        physics->se3                    = Se3r(Vector3r(0,0,0),q);
        physics->beamGroupMask          = beamGroupMask;
        physics->nodeGroupMask          = nodeGroupMask;
        physics->useBendTensileSoftening= useBendTensileSoftening;
        physics->useStiffnessSoftening  = useStiffnessSoftening;
        
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
	
 
void LatticeExample::imposeTranslation(shared_ptr<Scene>& rootBody, Vector3r min, Vector3r max, Vector3r direction, Real displacement,Vector3r blocked_xyz)
{
	shared_ptr<DisplacementEngine> translationCondition = shared_ptr<DisplacementEngine>(new DisplacementEngine);
 	translationCondition->displacement  = displacement;
	direction.Normalize();
 	translationCondition->translationAxis = direction;

        // FIXME: WTF ???
        rootBody->engines.push_back((rootBody->engines)[rootBody->engines.size()-1]);
        (rootBody->engines)[rootBody->engines.size()-2]=(rootBody->engines)[rootBody->engines.size()-3];
        (rootBody->engines)[rootBody->engines.size()-3]=(rootBody->engines)[rootBody->engines.size()-4];
        (rootBody->engines)[rootBody->engines.size()-4]=(rootBody->engines)[rootBody->engines.size()-5];
        (rootBody->engines)[rootBody->engines.size()-5]=(rootBody->engines)[rootBody->engines.size()-6];
        (rootBody->engines)[rootBody->engines.size()-6]=(rootBody->engines)[rootBody->engines.size()-7];
        (rootBody->engines)[rootBody->engines.size()-7]=translationCondition;
        translationCondition->subscribedBodies.clear();
        
        BodyContainer::iterator bi    = rootBody->bodies->begin();
	BodyContainer::iterator biEnd = rootBody->bodies->end();

	if(region_single_node_ABCDEF)
	{
		Vector3r MIN(min);
		Vector3r MAX(max);
		Vector3r center = (MIN+MAX)*0.5;
		int best=-1;
		Real dist = (MIN-MAX).Length();
		Real prev_dist = dist;
		for(int run=0 ; run < 3 ; ++run )
		{
			center = (MIN+MAX)*0.5;
			dist = (MIN-MAX).Length();
			prev_dist = dist;
			bi    = rootBody->bodies->begin();
			biEnd = rootBody->bodies->end();
			for(  ; bi!=biEnd ; ++bi )
			{
				shared_ptr<Body> b = *bi;
			
				if( b->getGroupMask() & nodeGroupMask )
				{
					Vector3r pos = b->physicalParameters->se3.position;
					if(        pos[0] >= min[0] 
						&& pos[1] >= min[1] 
						&& pos[2] >= min[2] 
						&& pos[0] <= max[0] 
						&& pos[1] <= max[1] 
						&& pos[2] <= max[2] 
						&& (b->getGroupMask() & nodeGroupMask)
						)
					{
					//	b->isDynamic = false;
						b->geometricalModel->diffuseColor = Vector3r(2.0,2.0,0.0);
					//	translationCondition->subscribedBodies.push_back(b->getId());

						dist = (pos-center).Length();
						if(dist < prev_dist)
						{
							best = b->getId();
							prev_dist = dist;
						}
						if(run == 0)
						{
							MIN=pos;
							MAX=pos;
							run = 1;
						}
						if(run == 1)
						{
							MIN=componentMinVector(MIN,pos);
							MAX=componentMaxVector(MAX,pos);
						}
					}
				}
			}
		}
		if(best != -1)
		{
			std::cerr << "INFO: single node in region, best ID is: " << best << ", dist: " << prev_dist 
				<< ", at coords: " << (*(rootBody->bodies))[best]->physicalParameters->se3.position << "\n";
			translationCondition->subscribedBodies.push_back(best);
			(*(rootBody->bodies))[best]->geometricalModel->diffuseColor = Vector3r(0.0,0.0,3.0);
			(*(rootBody->bodies))[best]->isDynamic = false;
			(*(rootBody->bodies))[best]->physicalParameters->setDOFfromVector3r(blocked_xyz);
		}
		else
		{
			std::cerr << "WARNING: cannot find a single node in this region!\n";
		}
	}
	else
	{
		for(  ; bi!=biEnd ; ++bi )
		{
			shared_ptr<Body> b = *bi;
		
			if( b->getGroupMask() & nodeGroupMask )
			{
				Vector3r pos = b->physicalParameters->se3.position;
				if(        pos[0] >= min[0] 
					&& pos[1] >= min[1] 
					&& pos[2] >= min[2] 
					&& pos[0] <= max[0] 
					&& pos[1] <= max[1] 
					&& pos[2] <= max[2] 
					&& (b->getGroupMask() & nodeGroupMask)
					)
				{
					b->isDynamic = false;
					b->geometricalModel->diffuseColor = Vector3r(2.0,2.0,0.0);
					translationCondition->subscribedBodies.push_back(b->getId());
					b->physicalParameters->setDOFfromVector3r(blocked_xyz);
				}
			}
		}
	}
}

bool LatticeExample::isDeleted(Vector3r pos, Vector3r min, Vector3r max)
{
	return(
		   pos[0] > min[0] 
		&& pos[1] > min[1] 
		&& pos[2] > min[2] 
		&& pos[0] < max[0] 
		&& pos[1] < max[1] 
		&& pos[2] < max[2] 
	      );
}

Real sectionArea(Vector3r min,Vector3r max,Vector3r box)
{
	Vector3r mi = componentMaxVector(Vector3r(0,0,0),min);
	Vector3r Ma = componentMinVector(box            ,max);
	if(mi[0]<Ma[0] && mi[1]<Ma[1])
		return (Ma[0]-mi[0])*(Ma[1]-mi[1]);
	else
		return 0;
}

Real sectionVolu(Vector3r min,Vector3r max,Vector3r box)
{
	Vector3r mi = componentMaxVector(Vector3r(0,0,0),min);
	Vector3r Ma = componentMinVector(box            ,max);
	if(mi[0]<Ma[0] && mi[1]<Ma[1] && mi[2]<Ma[2])
		return (Ma[0]-mi[0])*(Ma[1]-mi[1])*(Ma[2]-mi[2]);
	else
		return 0;
}

Real LatticeExample::deletedArea()
{
	return
		  sectionArea(regionDelete_A_min,regionDelete_A_max,speciemen_size_in_meters)
        	+ sectionArea(regionDelete_B_min,regionDelete_B_max,speciemen_size_in_meters)
        	+ sectionArea(regionDelete_C_min,regionDelete_C_max,speciemen_size_in_meters)
        	+ sectionArea(regionDelete_D_min,regionDelete_D_max,speciemen_size_in_meters)
        	+ sectionArea(regionDelete_E_min,regionDelete_E_max,speciemen_size_in_meters)
        	+ sectionArea(regionDelete_F_min,regionDelete_F_max,speciemen_size_in_meters)

		+ sectionArea(regionDelete_1_min,regionDelete_1_max,speciemen_size_in_meters)
        	+ sectionArea(regionDelete_2_min,regionDelete_2_max,speciemen_size_in_meters)
        	+ sectionArea(regionDelete_3_min,regionDelete_3_max,speciemen_size_in_meters)
        	+ sectionArea(regionDelete_4_min,regionDelete_4_max,speciemen_size_in_meters)
        	+ sectionArea(regionDelete_5_min,regionDelete_5_max,speciemen_size_in_meters);
};

Real LatticeExample::deletedVolume()
{
	return
		  sectionVolu(regionDelete_A_min,regionDelete_A_max,speciemen_size_in_meters)
        	+ sectionVolu(regionDelete_B_min,regionDelete_B_max,speciemen_size_in_meters)
        	+ sectionVolu(regionDelete_C_min,regionDelete_C_max,speciemen_size_in_meters)
        	+ sectionVolu(regionDelete_D_min,regionDelete_D_max,speciemen_size_in_meters)
        	+ sectionVolu(regionDelete_E_min,regionDelete_E_max,speciemen_size_in_meters)
        	+ sectionVolu(regionDelete_F_min,regionDelete_F_max,speciemen_size_in_meters)

		+ sectionVolu(regionDelete_1_min,regionDelete_1_max,speciemen_size_in_meters)
        	+ sectionVolu(regionDelete_2_min,regionDelete_2_max,speciemen_size_in_meters)
        	+ sectionVolu(regionDelete_3_min,regionDelete_3_max,speciemen_size_in_meters)
        	+ sectionVolu(regionDelete_4_min,regionDelete_4_max,speciemen_size_in_meters)
        	+ sectionVolu(regionDelete_5_min,regionDelete_5_max,speciemen_size_in_meters);
};

bool LatticeExample::notDeleted(Vector3r pos)
{
//	return true;
	return (!(
		   isDeleted(pos,regionDelete_A_min,regionDelete_A_max)
        	|| isDeleted(pos,regionDelete_B_min,regionDelete_B_max)
        	|| isDeleted(pos,regionDelete_C_min,regionDelete_C_max)
        	|| isDeleted(pos,regionDelete_D_min,regionDelete_D_max)
        	|| isDeleted(pos,regionDelete_E_min,regionDelete_E_max)
        	|| isDeleted(pos,regionDelete_F_min,regionDelete_F_max)

		|| isDeleted(pos,regionDelete_1_min,regionDelete_1_max)
        	|| isDeleted(pos,regionDelete_2_min,regionDelete_2_max)
        	|| isDeleted(pos,regionDelete_3_min,regionDelete_3_max)
        	|| isDeleted(pos,regionDelete_4_min,regionDelete_4_max)
        	|| isDeleted(pos,regionDelete_5_min,regionDelete_5_max)
	));
};

void LatticeExample::regionDelete(shared_ptr<Scene>& rootBody, Vector3r min, Vector3r max)
{
        vector<unsigned int> futureDeletes;
        
        BodyContainer::iterator bi    = rootBody->bodies->begin();
        BodyContainer::iterator biEnd = rootBody->bodies->end();
        for(  ; bi!=biEnd ; ++bi )
        {
                shared_ptr<Body> b = *bi;
        
                if( b->getGroupMask() & beamGroupMask )
                {
                        Vector3r pos = b->physicalParameters->se3.position;
                        if(        pos[0] > min[0] 
                                && pos[1] > min[1] 
                                && pos[2] > min[2] 
                                && pos[0] < max[0] 
                                && pos[1] < max[1] 
                                && pos[2] < max[2] 
                                )
                                futureDeletes.push_back( b->getId() );
                }
        }
        
        vector<unsigned int>::iterator vend = futureDeletes.end();
        for( vector<unsigned int>::iterator vsta = futureDeletes.begin() ; vsta != vend ; ++vsta)
                rootBody->bodies->erase(*vsta); 
}

void LatticeExample::nonDestroy(shared_ptr<Scene>& rootBody, Vector3r min, Vector3r max)
{
	std::list<unsigned int> marked;
        
        BodyContainer::iterator bi    = rootBody->bodies->begin();
        BodyContainer::iterator biEnd = rootBody->bodies->end();
        for(  ; bi!=biEnd ; ++bi )
        {
                shared_ptr<Body> b = *bi;
        
                if( b->getGroupMask() & beamGroupMask )
                {
                        Vector3r pos = b->physicalParameters->se3.position;
                        if(        pos[0] > min[0] 
                                && pos[1] > min[1] 
                                && pos[2] > min[2] 
                                && pos[0] < max[0] 
                                && pos[1] < max[1] 
                                && pos[2] < max[2] 
                                )
                                marked.push_back( b->getId() );
                }
        }
        
	std::list<unsigned int>::iterator vend = marked.end();
        for( std::list<unsigned int>::iterator vsta = marked.begin() ; vsta != vend ; ++vsta)
        {
                LatticeBeamParameters* beam = static_cast<LatticeBeamParameters*>( ((*(rootBody->bodies))[*vsta])->physicalParameters.get());
                beam->criticalTensileStrain     = 0.9;
                beam->criticalCompressiveStrain = 0.9;
                beam->longitudalStiffness       = nonDestroy_stiffness;
                beam->bendingStiffness          = nonDestroy_stiffness*0.1;
                beam->torsionalStiffness        = nonDestroy_stiffness*0.1;
//                (*(rootBody->bodies))[beam->id1]->geometricalModel->diffuseColor = Vector3r(0.2,0.5,0.7);
//                (*(rootBody->bodies))[beam->id2]->geometricalModel->diffuseColor = Vector3r(0.2,0.5,0.7);
        }
}


void LatticeExample::modifyCT(shared_ptr<Scene>& rootBody, Vector3r min, Vector3r max)
{
	std::list<unsigned int> marked;
        
        BodyContainer::iterator bi    = rootBody->bodies->begin();
        BodyContainer::iterator biEnd = rootBody->bodies->end();
        for(  ; bi!=biEnd ; ++bi )
        {
                shared_ptr<Body> b = *bi;
        
                if( b->getGroupMask() & beamGroupMask )
                {
                        Vector3r pos = b->physicalParameters->se3.position;
                        if(        pos[0] > min[0] 
                                && pos[1] > min[1] 
                                && pos[2] > min[2] 
                                && pos[0] < max[0] 
                                && pos[1] < max[1] 
                                && pos[2] < max[2] 
                                )
                                marked.push_back( b->getId() );
                }
        }
        
	std::list<unsigned int>::iterator vend = marked.end();
        for( std::list<unsigned int>::iterator vsta = marked.begin() ; vsta != vend ; ++vsta)
        {
                LatticeBeamParameters* beam = static_cast<LatticeBeamParameters*>( ((*(rootBody->bodies))[*vsta])->physicalParameters.get());
                beam->criticalTensileStrain     = CT;
        }
}

/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

int LatticeExample::overlapCount(std::vector<Circle>& c,float threshold)
{
	int res=0;
	for(size_t i=0 ; i<c.size() ; ++i)
		if(overlaps(c[i],c,true,threshold))
			++res;
	return res;
};

bool LatticeExample::overlaps(Circle& cc,std::vector<Circle>& c,bool no_self_overlap,float threshold)
{
	// check with circles
        std::vector<Circle>::iterator end=c.end();
        for(std::vector<Circle>::iterator i=c.begin();i!=end;++i)
        {
                float dist2 = std::pow(i->x - cc.x ,2)+std::pow(i->y - cc.y,2)+std::pow(i->z - cc.z,2);

		if(no_self_overlap && (dist2 == 0) && ((i->x==cc.x) && (i->y==cc.y) && (i->z==cc.z) && (i->d==cc.d)))
			continue;

                float r2    = std::pow( threshold*(i->d+cc.d)/2.0 ,2); // FIXME - 1.1 is hardcoded. van Mier's min distance is 1.1*(D1+D2)/2
                if(dist2<r2)
                        return true;
        }

	// check with fibres
	if(fibre_count > 0)
	{
		for(int i = 0 ; i < fibre_count ; ++i)
		{
			Vector3r pos = fibres[i].first;
			Vector3r del = fibres[i].second;
			for(int j = 0 ; j < beams_per_fibre ; ++j)
			{
				Vector3r p(
					pos[0] + 1.0*j*del[0],
					pos[1] + 1.0*j*del[1],
					pos[2] + 1.0*j*del[2] );
				Vector3r dist;
				Vector3r a(cc.x,cc.y,cc.z);
				dist = p-a;
				if(dist.Length() < cc.d*0.5 )
					return true;
			}
		}
	}

	// check if it's outside the concrete
//	if(no_Agg_outside)
//	{
//		Real AGGREGATES_X=speciemen_size_in_meters[0];
//	        Real AGGREGATES_Y=speciemen_size_in_meters[1];
//	        Real AGGREGATES_Z=speciemen_size_in_meters[2];
//		for(Real A=std::max(0.0,cc.x-cc.d*0.5) ; A<=std::min(AGGREGATES_X,cc.x+cc.d*0.5) ; A+=cellsizeUnit_in_meters*0.2)
//		for(Real B=std::max(0.0,cc.y-cc.d*0.5) ; B<=std::min(AGGREGATES_Y,cc.y+cc.d*0.5) ; B+=cellsizeUnit_in_meters*0.2)
//		for(Real C=std::max(0.0,cc.z-cc.d*0.5) ; C<=std::min(AGGREGATES_Z,cc.z+cc.d*0.5) ; C+=cellsizeUnit_in_meters*0.2)
//			if(! notDeleted(Vector3r(A,B,C)) )
//			{
//				if( std::pow(cc.x - A,2.0) + std::pow(cc.y - B,2.0) + std::pow(cc.z - C,2.0) < cc.d*cc.d*0.25 )
//					return true;
//			}
//	}

        return false;
};

int LatticeExample::aggInside(Vector3r& a,Vector3r& b,std::vector<Circle>& c, Real cellsizeUnit_in_meters)
{ // checks if nodes 'a','b' are inside any of aggregates from list 'c'
        int res=0;
        std::vector<Circle>::iterator end=c.end();
        for(std::vector<Circle>::iterator i=c.begin();i!=end;++i)
        {
        //      if(i->r < cellsizeUnit_in_meters) // FIXME
        //              continue;

                float dist2 = std::pow(i->x - a[0],2)+std::pow(i->y - a[1],2)+std::pow(i->z - a[2],2);
                float r2    = std::pow(i->d*0.5,2);
                if(dist2<r2) res=1; else res=0;

                dist2 = std::pow(i->x - b[0],2)+std::pow(i->y - b[1],2)+std::pow(i->z - b[2],2);
                if(dist2<r2) ++res;

                if(res!=0) return res;
        }
        return false;
}

float LatticeExample::aggsAreas(std::vector<Circle>& c)
{
        float aggArea=0.0;
        std::vector<Circle>::iterator end=c.end();
        for(std::vector<Circle>::iterator i=c.begin();i!=end;++i)
                aggArea += 3.14159265358979323846*std::pow(i->d*0.5 ,2);
        return aggArea;
}

float LatticeExample::aggsVolumes(std::vector<Circle>& c)
{
        float aggVolume=0.0;
        std::vector<Circle>::iterator end=c.end();
        for(std::vector<Circle>::iterator i=c.begin();i!=end;++i)
                aggVolume += (4.0/3.0)*3.14159265358979323846*std::pow(i->d*0.5 ,3);
        return aggVolume;
}

// random
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>

void LatticeExample::addAggregates(shared_ptr<Scene>& rootBody)
{
        // first make a list of circles
        std::vector<Circle> c;
        std::vector<int> penalty;

        Real AGGREGATES_X=speciemen_size_in_meters[0];
        Real AGGREGATES_Y=speciemen_size_in_meters[1];
        Real AGGREGATES_Z=speciemen_size_in_meters[2];
	if(AGGREGATES_Z < cellsizeUnit_in_meters )
		AGGREGATES_Z = 0.0;

        typedef boost::minstd_rand StdGenerator;
        static StdGenerator generator;
        static boost::variate_generator<StdGenerator&, boost::uniform_real<> >
                random1(generator, boost::uniform_real<>(0,1));
///////////////////////////////////////////////////////////////
/////// stage 2 ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////

        Real MAX_DIAMETER  =aggregateMaxDiameter;
        Real MIN_DIAMETER  =aggregateMinDiameter;
        Real MEAN_DIAMETER =aggregateMeanDiameter;
        Real SIGMA_DIAMETER=aggregateSigmaDiameter;
        static boost::variate_generator<StdGenerator&, boost::normal_distribution<> > 
                randomN(generator, boost::normal_distribution<>(MEAN_DIAMETER,SIGMA_DIAMETER));

        std::cerr << "generating aggregates ... ";
	setStatus(  "generating aggregates...");
        do
        {
		if(shouldTerminate())
		{
			setTerminate(false);
			break;
		}

                Circle cc;
                cc.x=random1()*AGGREGATES_X, cc.y=random1()*AGGREGATES_Y, cc.z=AGGREGATES_Z==0?0:random1()*AGGREGATES_Z;
                do { cc.d=randomN(); } while (cc.d>=MAX_DIAMETER || cc.d<=MIN_DIAMETER);
                for(int i=0 ; i<1000 ; ++i)
		/* IGNORE OVERLAPPING FOR NOW, just add as much as we need to!
                        if(overlaps(cc,c))
                                cc.x=random1()*AGGREGATES_X, cc.y=random1()*AGGREGATES_Y, cc.z=AGGREGATES_Z==0?0:random1()*AGGREGATES_Z;
                        else
		*/
                        {
                                c.push_back(cc);
                //              std::cerr << cc.x << " " << cc.y << " " << cc.d << "\n";
                                break;
                        }
		if(AGGREGATES_Z == 0)
			setProgress((aggsAreas(c)/(AGGREGATES_X*AGGREGATES_Y - ((no_Agg_outside)?(deletedArea()):(0.0)) ))/(aggregatePercent/100.0));
		else
			setProgress((aggsVolumes(c)/(AGGREGATES_X*AGGREGATES_Y*AGGREGATES_Z - ((no_Agg_outside)?(deletedVolume()):(0.0)) ))/(aggregatePercent/100.0));
        }
        //while(aggregatePercent/100.0 > aggsAreas(c)/(AGGREGATES_X*AGGREGATES_Y) );
        while( progress() < 1.0 );

	std::cerr << "placing aggregates with repulsion ... ";
	setStatus(   "aggregates repulsion ...");
	setProgress(0);
	penalty.clear();
	penalty.resize(c.size(),0);

///////////////////////////////////////////////////////////////
/////// stage 3 ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////
        
	
	int overlap_count;
	overlap_count=overlapCount(c,1.05);
	int begin_overlap_count(overlap_count);

// repulsion !!
	for( ; overlap_count > 0 ; )
	{
		overlap_count=overlapCount(c,1.05);

		int last_overlap( ((int)(overlap_count)/10)*10+10 );
		if(overlap_count <= last_overlap - 10)
		{
			//setStatus(std::string("repulsion ")+boost::lexical_cast<std::string>(overlap_count));
			last_overlap = ((int)(overlap_count)/10)*10;
			std::cerr << " . " << last_overlap;
			if(last_overlap <= 11)
				last_overlap = 11;
		}
		setProgress(1.0 - (float)(overlap_count)/(float)(begin_overlap_count));

		std::vector<Vector3r > moves;
		moves.clear();
		for(unsigned int i = 0 ; i < c.size() ; ++i )
		{
			Vector3r d(0,0,0);

			Vector3r c1(c[i].x,c[i].y,c[i].z);

			//emulate periodic boundary
			//for(int px = -1 ; px < 2 ; ++px )
			//for(int py = -1 ; py < 2 ; ++py )
			//for(int pz = ((AGGREGATES_Z==0)?0:-1) ; pz < ((AGGREGATES_Z==0)?1:2) ; ++pz )
			int px(0),py(0),pz(0);
			{
				Vector3r PERIODIC_DELTA(px*AGGREGATES_X,py*AGGREGATES_Y,pz*AGGREGATES_Z);
				for(unsigned int j = 0 ; j < c.size() ; ++j )
					if(i != j)
					{
						Vector3r c2 = Vector3r(c[j].x, c[j].y, c[j].z) + PERIODIC_DELTA;

						Vector3r dir=c1-c2;
						Real r = dir.Normalize(); // dir is unit vector, r is a distance
						if(r < cellsizeUnit_in_meters)
						{
							r=cellsizeUnit_in_meters;
							dir[0]=random1()-0.5, dir[1]=random1()-0.5, dir[2]=((AGGREGATES_Z==0)?(0):(random1()-0.5));
							dir.Normalize();
						}
						// weak distance repulsion
						d += dir * 1/(r*r)*0.3;

						// strong repulsion when they are intersecting
						if(2*r < (c[i].d+c[j].d)*1.11)
						{
							d += (dir * 1/(r*r))*3;
							++penalty[i];
							++penalty[j];
						}
					}
			}

			// repulsion from walls.
			Vector3r MAX(AGGREGATES_X, AGGREGATES_Y, AGGREGATES_Z);
			for(int I=0 ; I<((AGGREGATES_Z==0)?2:3) ; ++I)
			{
				//if(c1[I] > 0 && c1[I] < MAX[I])
				//	d[I] += (1/(c1[I]*c1[I]) - 1/((MAX[I]-c1[I])*(MAX[I]-c1[I])))*0.5;

				if(c1[I] > 0 && c1[I] < c[i].d*0.5 )
					d[I] += (1/(c1[I]*c1[I]));
				else
				if(c1[I] < MAX[I] && c1[I] > MAX[I]-c[i].d*0.5)
					d[I] += (-1/((MAX[I]-c1[I])*(MAX[I]-c1[I])));
			}

			// check with fibres
			if(fibre_count > 0)
			{
				for(int I = 0 ; I < fibre_count ; ++I)
				{
					Vector3r pos = fibres[I].first;
					Vector3r del = fibres[I].second;
					for(int J = 0 ; J < beams_per_fibre ; ++J)
					{
						Vector3r c2(
							pos[0] + 1.0*J*del[0],
							pos[1] + 1.0*J*del[1],
							pos[2] + 1.0*J*del[2] );
						
						Vector3r dir=c1-c2;
						
						Real r = dir.Normalize(); // dir is unit vector, r is a distance
						if(r < cellsizeUnit_in_meters)
						{
							r=cellsizeUnit_in_meters;
							dir[0]=random1()-0.5, dir[1]=random1()-0.5, dir[2]=((AGGREGATES_Z==0)?(0):(random1()-0.5));
							dir.Normalize();
						}
						// weak repulnsion
						d += (dir * 1/(r*r))*(0.1/beams_per_fibre);

						// strong repulsion when they are intersecting
						if(2*r < (c[i].d)*1.11)
						{
							d += (dir * 1/(r*r));
							++penalty[i];
						}
					}
				}
			}

			moves.push_back(d);
		}

		assert(moves.size() == c.size() );
		Real maxl=0;
		for(unsigned int i = 0 ; i < moves.size() ; ++i )
			maxl = std::max(moves[i].Length(),maxl);
		for(unsigned int i = 0 ; i < moves.size() ; ++i )
			moves[i] = cellsizeUnit_in_meters*0.3*moves[i]/(maxl==0?1:maxl);
	
		for(unsigned int i = 0 ; i < moves.size() ; ++i )
		{
			c[std::floor(((float)(i)))].x+=moves[i][0];
			c[std::floor(((float)(i)))].y+=moves[i][1];
			c[std::floor(((float)(i)))].z+=moves[i][2];
		}
	
		for(unsigned int i = 0 ; i < c.size() ; ++i )
		{
			Vector3r c1(c[i].x, c[i].y, c[i].z);
			if(   c1[0] < 0 
			   || c1[1] < 0 
			   || c1[2] < 0
			   || c1[0] > AGGREGATES_X
			   || c1[1] > AGGREGATES_Y
			   || c1[2] > AGGREGATES_Z
			   || penalty[i] > (beams_per_fibre+10)*1500)
			{
				//std::cerr << "aggregate: putting again randomly ";
				float limit(1.31);
				int count(0);
				do
				{
					//std::cerr << ".";
					c[i].x = random1()*AGGREGATES_X;
					c[i].y = random1()*AGGREGATES_Y;
					c[i].z = ((AGGREGATES_Z==0)?(0):(random1()*AGGREGATES_Z));
					++count;
					if(count>1000 && limit > 0.85)
					{
						limit -= 0.02;
						count=0;
						//std::cerr << "\nlimit set to: " << limit << " ";
					}
				} while(overlaps(c[i],c,true,limit));
				penalty[i]=0;
				//std::cerr << "\n";
			}
		}
		
	//std::cerr << "frame= " << frame++ << "\n";

		if(shouldTerminate()) return;
	}

        std::cerr << "done. " << c.size() << " area: " << aggsAreas(c)/(AGGREGATES_X*AGGREGATES_Y) << " vol: " << aggsVolumes(c)/(AGGREGATES_X*AGGREGATES_Y*AGGREGATES_Z) << "\n";

/*
	if(no_Agg_outside)
	{
		for(Real A=0 ; A<=AGGREGATES_X ; A+=cellsizeUnit_in_meters*0.3)
		for(Real B=0 ; B<=AGGREGATES_Y ; B+=cellsizeUnit_in_meters*0.3)
		for(Real C=0 ; C<=AGGREGATES_Z ; C+=cellsizeUnit_in_meters*0.3)
			if(! notDeleted(Vector3r(A,B,C)) )
			{ // remove a circle if it has inside Vector3r(A,B,C)
				std::vector<int> c_flag;c_flag.clear();c_flag.resize(c.size(),1);
				for(int i=0;i<c.size();++i)
					if( std::pow(c[i].x - A,2.0) + std::pow(c[i].y - B,2.0) + std::pow(c[i].z - C,2.0) < c[i].d*c[i].d*0.25 )
						c_flag[i]=0;

				std::vector<Circle> c_good;c_good.clear();
				for(int i=0;i<c.size();++i)
					if(c_flag[i]==1)
						c_good.push_back(c[i]);
				c=c_good;
			}
	}
*/

        { // set different properties for beams that lie in an aggregate
          // parametrize from above - takes three arguments: 
          // - MaterialParameters of aggregate, 
          // - MaterialParameters of bond, 
          // - list of circles.
                BodyContainer::iterator bi    = rootBody->bodies->begin();
                BodyContainer::iterator biEnd = rootBody->bodies->end();
                //float all_bodies = rootBody->bodies->size();
                //int current = 0;
                for(  ; bi!=biEnd ; ++bi )  // loop over all beams
                {
                //      if( ++current % 100 == 0 )
                //              cerr << "aggregate beams: " << current << " , " << ((static_cast<float>(current)/all_bodies)*100.0) << " %\n";
                                
                        Body* body = (*bi).get();
                        if( ! ( body->getGroupMask() & beamGroupMask ) )
                                continue; // skip non-beams
			
                        LatticeBeamParameters* beam     = static_cast<LatticeBeamParameters*>(body->physicalParameters.get());
                        int ovv = aggInside(     (*(rootBody->bodies))[beam->id1]->physicalParameters->se3.position
                                                ,(*(rootBody->bodies))[beam->id2]->physicalParameters->se3.position
                                                ,c, cellsizeUnit_in_meters );
                        if(ovv==2) // aggregate
                        {
                                beam->longitudalStiffness       = agg_longStiffness_noUnit;
                                beam->bendingStiffness          = agg_bendStiffness_noUnit;
                                beam->torsionalStiffness        = agg_torsStiffness_noUnit;
                                beam->criticalTensileStrain     = agg_critTensileStrain;
                                beam->criticalCompressiveStrain = agg_critCompressStrain;
			
				body->groupMask = beamGroupMask + aggGroupMask;
                
//                                (*(rootBody->bodies))[beam->id1]->geometricalModel->diffuseColor = Vector3r(0.6,0.2,0.0);
//                                (*(rootBody->bodies))[beam->id2]->geometricalModel->diffuseColor = Vector3r(0.6,0.2,0.0);
                        }
                        else if(ovv==1) // bond
                        {
                                beam->longitudalStiffness       = bond_longStiffness_noUnit;
                                beam->bendingStiffness          = bond_bendStiffness_noUnit;
                                beam->torsionalStiffness	= bond_torsStiffness_noUnit;
                                beam->criticalTensileStrain     = bond_critTensileStrain;
                                beam->criticalCompressiveStrain = bond_critCompressStrain;
                
				body->groupMask = beamGroupMask + bondGroupMask;
                
//                                (*(rootBody->bodies))[beam->id1]->geometricalModel->diffuseColor = Vector3r(0.6,0.6,0.0);
//                                (*(rootBody->bodies))[beam->id2]->geometricalModel->diffuseColor = Vector3r(0.6,0.6,0.0);
                        }
			else // matrix
				matrix_total+=1.0;
                }
        }
}

//
// old method - just put randomly
//
//void LatticeExample::makeFibres()
//{
//	fibres.clear();
//
//        Real AGGREGATES_X=speciemen_size_in_meters[0];
//        Real AGGREGATES_Y=speciemen_size_in_meters[1];
//        Real AGGREGATES_Z=speciemen_size_in_meters[2];
//	if(AGGREGATES_Z < cellsizeUnit_in_meters )
//		AGGREGATES_Z = 0.0;
//
//        typedef boost::minstd_rand StdGenerator;
//        static StdGenerator generator;
//        static boost::variate_generator<StdGenerator&, boost::uniform_real<> >
//                random1(generator, boost::uniform_real<>(0,1));
//
//	for(int i = 0 ; i < fibre_count ; ++i)
//        {
//                Vector3r cc;
//		Vector3r del;
//                cc[0]=random1()*AGGREGATES_X, cc[1]=random1()*AGGREGATES_Y, cc[2]=AGGREGATES_Z==0?0:random1()*AGGREGATES_Z;
//                del[0]=random1()-0.5, del[1]=random1()-0.5, del[2]=((AGGREGATES_Z==0)?(0):(random1()-0.5));
//		del.Normalize();
//		del=cellsizeUnit_in_meters * del;
//		fibres.push_back(std::make_pair(cc,del));
//        }
//}


// new method - equally balance fibres over volume using repulsion
void LatticeExample::makeFibres()
{
        Real AGGREGATES_X=speciemen_size_in_meters[0];
        Real AGGREGATES_Y=speciemen_size_in_meters[1];
        Real AGGREGATES_Z=speciemen_size_in_meters[2];
	if(AGGREGATES_Z < cellsizeUnit_in_meters )
		AGGREGATES_Z = 0.0;

        typedef boost::minstd_rand StdGenerator;
        static StdGenerator generator;
        static boost::variate_generator<StdGenerator&, boost::uniform_real<> >
                random1(generator, boost::uniform_real<>(0,1));

///////////////////////////////////////////////////////////////
/////// stage 0 ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////

	setStatus("balancing fibres...");
	std::cerr << "fibres: ";
	fibres.clear();

	for(int i = 0 ; i < fibre_count ; ++i)
        {
                Vector3r cc;
		Vector3r del;
                cc[0]=random1()*AGGREGATES_X, cc[1]=random1()*AGGREGATES_Y, cc[2]=((AGGREGATES_Z==0)?(0):(random1()*AGGREGATES_Z));

                //del[0]=random1()-0.5, del[1]=random1()-0.5, del[2]=((AGGREGATES_Z==0)?(0):(random1()-0.5));

		if(fibres_horizontal || fibres_vertical)
		{
			if(fibres_horizontal)
			{
				std::cerr << "fibres_horizontal\n";

				del[0]=0;
				while(del[0]*del[0] < 0.15) del[0]=random1()-0.5;

				del[1]=(random1()-0.5)*0.2;

				del[2]=((AGGREGATES_Z==0)?(0):(random1()-0.5));
			}
			else
			{
				std::cerr << "fibres_vertical\n";

				del[0]=(random1()-0.5)*0.2;

				del[1]=0;
				while(del[1]*del[1] < 0.15) del[1]=random1()-0.5;

				del[2]=((AGGREGATES_Z==0)?(0):(random1()-0.5));
			}
		}
		else
		{
			del[0]=random1()-0.5;
			del[1]=random1()-0.5;
			del[2]=((AGGREGATES_Z==0)?(0):(random1()-0.5));
		}

		del.Normalize();
		del=cellsizeUnit_in_meters * del;
		fibres.push_back(std::make_pair(cc,del));
        }
/*
	   for(int i = 0 ; i < fibre_count ; ++i)
	   {
		   Vector3r cc;
		   Vector3r del;
		   cc[0]=random1()*AGGREGATES_X, cc[1]=random1()*AGGREGATES_Y, cc[2]=((AGGREGATES_Z==0)?(0):(random1()*AGGREGATES_Z));

		   del[0]=random1()-0.5;
		   del[1]=random1()-0.5;
		   del[2]=((AGGREGATES_Z==0)?(0):(random1()-0.5));

		   del.Normalize();
		   del=cellsizeUnit_in_meters * del;
		   fibres.push_back(std::make_pair(cc,del));
	   }
*/
///////////////////////////////////////////////////////////////
/////// stage 1 ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////

// repulsion !!
	for(int frame=0; frame < fibre_balancing_iterations ; ++frame)
	{
		std::vector<Vector3r > moves;
		moves.clear();
		for(unsigned int i = 0 ; i < fibres.size() ; ++i )
		{
			Vector3r d(0,0,0);

			////Vector3r c1 = fibres[i].first + fibres[i].second*beams_per_fibre*0.5;

		for(float PART_1 = 0 ; PART_1 <= 1.0 ; PART_1 += 0.2 )
		{
			Vector3r c1 = fibres[i].first + fibres[i].second*beams_per_fibre*PART_1;

			//emulate periodic boundary
			//for(int px = -1 ; px < 2 ; ++px )
			//for(int py = -1 ; py < 2 ; ++py )
			//for(int pz = ((AGGREGATES_Z==0)?0:-1) ; pz < ((AGGREGATES_Z==0)?1:2) ; ++pz )
			int px(0),py(0),pz(0);
			{
				Vector3r PERIODIC_DELTA(px*AGGREGATES_X,py*AGGREGATES_Y,pz*AGGREGATES_Z);
				for(unsigned int j = 0 ; j < fibres.size() ; ++j )
					if(i != j)
					{
						////Vector3r c2 = fibres[j].first + fibres[j].second*beams_per_fibre*0.5 + PERIODIC_DELTA;


					for(float PART_2 = 0 ; PART_2 <= 1.0 ; PART_2 += 0.2 )
					{
						Vector3r c2 = fibres[j].first + fibres[j].second*beams_per_fibre*PART_2 + PERIODIC_DELTA;


						Vector3r dir=c1-c2;
						Real r = dir.Normalize(); // dir is unit vector, r is a distance
						if(r < cellsizeUnit_in_meters)
						{
							r=cellsizeUnit_in_meters;
							dir[0]=random1()-0.5, dir[1]=random1()-0.5, dir[2]=((AGGREGATES_Z==0)?(0):(random1()-0.5));
							dir.Normalize();
						}
						d += dir * 1/(r*r);

						// strong repulsion when they are intersecting (2D only)
						if(px==0 && py==0 && AGGREGATES_Z==0 && r<beams_per_fibre*cellsizeUnit_in_meters)
						{//stupid brute-force method
							bool overlaps=false;
							Vector3r start_i(fibres[i].first);
							Vector3r delta_i(fibres[i].second);

							Vector3r start_j(fibres[j].first);
							Vector3r delta_j(fibres[j].second);

							for(int I=0 ; I<beams_per_fibre ; ++I)
								for(int J=0 ; J<beams_per_fibre ; ++J)
								{
									float dist = ((start_i+delta_i*I) - (start_j+delta_j*J)).Length();
									if(dist < cellsizeUnit_in_meters*10.0)
										overlaps=true;
								}
							if(overlaps)
								d += (dir * 1/(r*r))*10;
						}
					}

					}
			}
			// repulsion from walls.
			Vector3r MAX(AGGREGATES_X, AGGREGATES_Y, AGGREGATES_Z);
			for(int I=0 ; I<((AGGREGATES_Z==0)?2:3) ; ++I)
			{
				//if(c1[I] > 0 && c1[I] < MAX[I])
				//	d[I] += (1/(c1[I]*c1[I]) - 1/((MAX[I]-c1[I])*(MAX[I]-c1[I])))*0.5;
				
				if(c1[I] > 0 && c1[I] < beams_per_fibre*cellsizeUnit_in_meters*0.2 )
					d[I] += (1/(c1[I]*c1[I]));
				else
				if(c1[I] < MAX[I] && c1[I] > MAX[I]-beams_per_fibre*cellsizeUnit_in_meters*0.2)
					d[I] += (-1/((MAX[I]-c1[I])*(MAX[I]-c1[I])));
			}
			moves.push_back(d);
		}

		}

		assert(moves.size() == fibres.size()*6 );
		Real maxl=0;
		for(unsigned int i = 0 ; i < moves.size() ; ++i )
			maxl = std::max(moves[i].Length(),maxl);
	//	std::cerr << "maxl= " << maxl << "\n";
		for(unsigned int i = 0 ; i < moves.size() ; ++i )
			moves[i] = cellsizeUnit_in_meters*moves[i]/maxl;
	
//		for(unsigned int i = 0 ; i < fibres.size() ; ++i )
//			fibres[i].first+=moves[i];
		for(unsigned int i = 0 ; i < moves.size() ; ++i )
			fibres[std::floor(((float)(i))/6.0f)].first+=moves[i];

	
		for(unsigned int i = 0 ; i < fibres.size() ; ++i )
		{
			Vector3r c1 = fibres[i].first + fibres[i].second*beams_per_fibre*0.5;
			if(   c1[0] < 0 
			   || c1[1] < 0 
			   || c1[2] < 0
			   || c1[0] > AGGREGATES_X
			   || c1[1] > AGGREGATES_Y
			   || c1[2] > AGGREGATES_Z)
			{
				//std::cerr << "putting again randomly\n";
				fibres[i].first = Vector3r(random1()*AGGREGATES_X, random1()*AGGREGATES_Y, ((AGGREGATES_Z==0)?(0):(random1()*AGGREGATES_Z)));
			}
		}
		
	//std::cerr << "frame= " << frame++ << "\n";

		if(shouldTerminate()) return;
		setProgress(1.0*frame/(1.0*fibre_balancing_iterations));
		int last_progress(0);
		if(last_progress == 0)
			std::cerr << "0% ";
		if( progress()*10 > last_progress )
			std::cerr << "... " << ++last_progress << "0% ";

	}
	std::cerr << "\n";

/*
	for(unsigned int i = 0 ; i < fibres.size() ; ++i )
		fibres[i].first+=Vector3r(random1()-0.5, random1()-0.5, ((AGGREGATES_Z==0)?(0):(random1()-0.5)) )*fibre_irregularity_noUnit*cellsizeUnit_in_meters;
*/
}


int LatticeExample::isFibre(Vector3r a,Vector3r b)
{
	int A=-1;
	int B=-1;
	for(int i = 0 ; i < fibre_count ; ++i)
        {
		Vector3r pos = fibres[i].first;
		Vector3r del = fibres[i].second;
		for(int j = 0 ; j < beams_per_fibre ; ++j)
		{
			Vector3r p(
				pos[0] + 1.0*j*del[0],
				pos[1] + 1.0*j*del[1],
				pos[2] + 1.0*j*del[2] );
			if(p == a) A=i;
			if(p == b) B=i;
		}
	}
	if(A+B==-2)
		return 0;
	if(A==B)
		return 2;
	return 1;
}

bool LatticeExample::fibreAllows(Vector3r a)
{
	for(int i = 0 ; i < fibre_count ; ++i)
        {
		Vector3r pos = fibres[i].first;
		Vector3r del = fibres[i].second / 3.0;
		for(int j = 0 ; j < beams_per_fibre * 3 ; ++j)
		{
			Vector3r p(
				pos[0] + 1.0*j*del[0],
				pos[1] + 1.0*j*del[1],
				pos[2] + 1.0*j*del[2] );
			Vector3r dist;
			dist = p-a;
			if(dist.Length() < cellsizeUnit_in_meters*fibre_allows )
				return false;
		}
	}
	return true;
}



void LatticeExample::makeFibreBeams(shared_ptr<Scene>& rootBody)
{
        { // set different properties for beams that are fibre
                BodyContainer::iterator bi    = rootBody->bodies->begin();
                BodyContainer::iterator biEnd = rootBody->bodies->end();
                for(  ; bi!=biEnd ; ++bi )  // loop over all beams
                {
                        Body* body = (*bi).get();
                        if( ! ( body->getGroupMask() & beamGroupMask ) )
                                continue; // skip non-beams

                        LatticeBeamParameters* beam     = static_cast<LatticeBeamParameters*>(body->physicalParameters.get());
			int fibreNodes=isFibre(   (*(rootBody->bodies))[beam->id1]->physicalParameters->se3.position
						 ,(*(rootBody->bodies))[beam->id2]->physicalParameters->se3.position);
			if(fibreNodes==2) // fibre
			{
				beam->longitudalStiffness       = fibre_longStiffness_noUnit;
                        	beam->bendingStiffness          = fibre_bendStiffness_noUnit;
                        	beam->torsionalStiffness	= fibre_torsStiffness_noUnit;
                        	beam->criticalTensileStrain     = fibre_critTensileStrain;
                        	beam->criticalCompressiveStrain = fibre_critCompressStrain;
                
//                        	(*(rootBody->bodies))[beam->id1]->geometricalModel->diffuseColor = Vector3r(2.0,2.0,0.0);
//                        	(*(rootBody->bodies))[beam->id2]->geometricalModel->diffuseColor = Vector3r(2.0,2.0,0.0);
				
				body->groupMask = beamGroupMask + fibreGroupMask;

				fibres_total+=1.0;
			}
			if(fibreNodes==1) // bond
			{
                                beam->longitudalStiffness       = fibre_bond_longStiffness_noUnit;
                                beam->bendingStiffness          = fibre_bond_bendStiffness_noUnit;
                                beam->torsionalStiffness	= fibre_bond_torsStiffness_noUnit;
                                beam->criticalTensileStrain     = fibre_bond_critTensileStrain;
                                beam->criticalCompressiveStrain = fibre_bond_critCompressStrain;
				
				body->groupMask = beamGroupMask + fibreBondGroupMask;

			}
                }
        }
}


YADE_PLUGIN((LatticeExample));

YADE_REQUIRE_FEATURE(PHYSPAR);

