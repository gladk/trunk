/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/core/FileGenerator.hpp>
#include<Wm3Vector3.h>
#include<yade/lib-base/yadeWm3.hpp>

class LatticeBeamParameters;
class StrainRecorder;

class LatticeExampleSimple : public FileGenerator
{
	private :
		int 		nodeGroupMask,beamGroupMask;

/////////////////////////////////////////////////////////////////////////
/////////////////////////  PARAMETERS ///////////////////////////////////
/////////////////////////////////////////////////////////////////////////
	// mesh generation
		Vector3r 	 speciemen_size_in_meters       // size
				,disorder_in_cellsize_unit;     // s
                Real             cell_size_in_meters            // g
				,max_rod_length_in_cellsize_unit

	// MaterialParameters of cement matrix
                                ,crit_tensile_strain            // E_min
                                ,crit_compr_strain              // E_max
                                ,longitudal_stiffness_no_unit   // k_l
                                ,bending_stiffness_no_unit      // k_b
                                ,torsional_stiffness_no_unit    // k_t
                                ;
	// conditions - displacement of top & bottom
		bool		 smooth_edges_displacement;
		Real		 y_top_bottom_displ_in_meters;

/////////////////////////////////////////////////////////////////////////
/////////////////////////  PARAMETERS END  //////////////////////////////
/////////////////////////////////////////////////////////////////////////

        public : 
		LatticeExampleSimple();
		virtual ~LatticeExampleSimple();

		bool generate();
		void make_simulation_loop(shared_ptr<Scene>& );
		void positionRootBody(shared_ptr<Scene>& rootBody);
		bool createNodeXYZ(shared_ptr<Body>& body, Vector3r position);
		void createBeam(shared_ptr<Body>& body, int i, int j);
		Real calcBeamPositionOrientationLength(shared_ptr<Body>& body, Real E_min, Real E_max, Real k_l, Real k_b, Real k_t);
		void create_angular_springs(int nodes_count,shared_ptr<Scene> rootBody);
		void calcBeamAngles(Body* body, BodyContainer* bodies, InteractionContainer* ints,std::vector< std::vector< unsigned int > >& connections);
		void calcAxisAngle(LatticeBeamParameters* beam1, BodyContainer* bodies, int otherId, InteractionContainer* ints, int thisId);
		void imposeTranslation(shared_ptr<Scene>& rootBody, Vector3r min, Vector3r max, Vector3r direction, Real displacement);

                REGISTER_CLASS_NAME(LatticeExampleSimple);
		REGISTER_BASE_CLASS_NAME(FileGenerator);
	REGISTER_ATTRIBUTES(FileGenerator,
		// mesh generation
		(speciemen_size_in_meters)     // size
		(cell_size_in_meters)          // g [m]        - cell size
		(disorder_in_cellsize_unit)    // s [-]        - disorder 

		// we can hide this parameter, because it's less important, just comment this line
		(max_rod_length_in_cellsize_unit)

		// MaterialParameters of cement matrix
		(crit_tensile_strain)           // E_min [%]    - default 0.02 %
		(crit_compr_strain)             // E_max [%]    - default 0.2 %
		(longitudal_stiffness_no_unit) // k_l [-]      - default 1.0
		(bending_stiffness_no_unit)    // k_b [-]      - default 0.6
		(torsional_stiffness_no_unit)  // k_t [-]      - default 0.6

		// conditions - displacement of top & bottom
		(smooth_edges_displacement)
		(y_top_bottom_displ_in_meters)
	)
};

REGISTER_SERIALIZABLE(LatticeExampleSimple);

