/*************************************************************************
*  Copyright (C) 2006 by Bruno Chareyre		                         *
*  bruno.chareyre@hmg.inpg.fr                                            *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/core/FileGenerator.hpp>

class ForceRecorder;
class AveragePositionRecorder;
class VelocityRecorder;
class TriaxialStressController;
class TriaxialCompressionEngine;
class TriaxialStateRecorder;
class GlobalStiffnessTimeStepper;

/*! \brief Isotropic compression + triaxial compression test

	TriaxialTest full documentation is here : http://yade-dem.org/wiki/TriaxialTest
	This preprocessor shows how to simulate a triaxial test in Yade. It is using the elastic-frictional contact law defined in ElasticContactLaw (similar procedures can be used for other contact laws).  It is designed to :
	1/ generate random loose packings and compress them under isotropic confining stress, either squeezing the packing between moving rigid boxes or expanding the particles while boxes are fixed (depending on flag "InternalCompaction").
	2/ simulate all sorts triaxial loading path (there is however a default loading path corresponding to constant lateral stress in 2 directions and constant strain rate on the third direction - this loading path is used when the flag AutoCompressionActivation = true, otherwise the simulation stops after isotropic compression)
	
	
	Essential engines :
	
	1/ The TrixialCompressionEngine is used for controlling the state of the sample and simulating loading paths. TrixialCompressionEngine inherits from TriaxialStressController, which can compute stress- strain-like quantities in the packing and maintain a constant level of stress at each boundary. TriaxialCompressionEngine has few more members in order to impose constant strain rate and control the transition between isotropic compression and triaxial test.
		
	2/ The class TriaxialStateRecorder is used to write to a file the history of stresses and strains.
	
	3/ TriaxialTest is currently using GlobalStiffnessTimeStepper to compute an appropriate dt for the numerical scheme. The TriaxialTest is the only preprocessor using these classes in Yade because they have been developed AFTER most of preprocessor examples, BUT they can be used in principle in any situation and they have nothing specifically related to the triaxial test.
	
	@note TriaxialStressController::ComputeUnbalancedForce(...) returns a value that can be useful for evaluating the stability of the packing. It is defined as (mean force on particles)/(mean contact force), so that it tends to 0 in a stable packing. This parameter is checked by TriaxialCompressionEngine to switch from one stage of the simulation to the next one (e.g. stop isotropic confinment and start axial loading)

	@note Compaction is done (1) by moving rigid boxes or (2) by increasing the sizes of the particles (decided using the option "internalCompaction" : true => size increase). Both algorithm needs numerical parameters to prevent instabilities. For instance, with method (1) maxWallVelocity is the maximum wall velocity, with method (2) finalMaxMultiplier is the max value of the multiplier applied on sizes at each iteration (always something like 1.00001). 
	
 */

class TriaxialTest : public FileGenerator
{
	private	:
		Vector3r	 gravity;
		Vector3r	 spheresColor;
		bool		 wall_top
				,wall_bottom
				,wall_1
				,wall_2
				,wall_3
				,wall_4
				,wall_top_wire
				,wall_bottom_wire
				,wall_1_wire
				,wall_2_wire
				,wall_3_wire
				,wall_4_wire
				,spheresRandomColor;
	
		shared_ptr<TriaxialCompressionEngine> triaxialcompressionEngine;
		shared_ptr<TriaxialStressController> triaxialstressController;
		shared_ptr<TriaxialStateRecorder> triaxialStateRecorder;
		shared_ptr<GlobalStiffnessTimeStepper> globalStiffnessTimeStepper;
			
		void createBox(shared_ptr<Body>& body, Vector3r position, Vector3r extents,bool wire);
		void createSphere(shared_ptr<Body>& body, Vector3r position, Real radius,bool big,bool dynamic);
		void createActors(shared_ptr<Scene>& scene);
		void positionRootBody(shared_ptr<Scene>& scene);
		typedef pair<Vector3r, Real> BasicSphere;	
	public : 
		~TriaxialTest ();
		bool generate(string& message);
		
		YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(
		TriaxialTest,FileGenerator,"Prepare a scene for triaxial tests. Full documentation in py/_extraDocs.py."
		,
   		((Vector3r,lowerCorner,Vector3r(0,0,0),,"Lower corner of the box."))
		((Vector3r,upperCorner,Vector3r(1,1,1),,"Upper corner of the box."))
		((string,importFilename,"",,"File with positions and sizes of spheres."))
		((string,Key,"",,"A code that is added to output filenames."))
		((string,fixedBoxDims,"",,"string that contains some subset (max. 2) of {'x','y','z'} ; contains axes will have box dimension hardcoded, even if box is scaled as mean_radius is prescribed: scaling will be applied on the rest."))
		((string,WallStressRecordFile,"./WallStresses"+Key,,""))					
		((bool,internalCompaction,false,,"flag for choosing between moving boundaries or increasing particles sizes during the compaction stage."))
		((bool,biaxial2dTest,false,,"FIXME : what is that?"))
		((bool,autoCompressionActivation,true,,"Do we just want to generate a stable packing under isotropic pressure (false) or do we want the triaxial loading to start automatically right after compaction stage (true)?"))
		((bool,autoUnload,true,,"auto adjust the isotropic stress state from :yref:`TriaxialTest::sigmaIsoCompaction` to :yref:`TriaxialTest::sigmaLateralConfinement` if they have different values. See docs for :yref:`TriaxialCompressionEngine::autoUnload`"))
		((bool,autoStopSimulation,false,,"freeze the simulation when conditions are reached (don't activate this if you want to be able to run/stop from Qt GUI)"))
		((bool,noFiles,false,,"Do not create any files during run (.xml, .spheres, wall stress records)"))
		((bool,facetWalls,false,,"Use facets for boundaries (not tested)"))
		((bool,wallWalls,false,,"Use walls for boundaries (not tested)"))
		
		((Real,thickness,0.001,,"thickness of boundaries. It is arbitrary and should have no effect"))
		((Real,maxMultiplier,1.01,,"max multiplier of diameters during internal compaction (initial fast increase)"))
		((Real,finalMaxMultiplier,1.001,,"max multiplier of diameters during internal compaction (secondary precise adjustment)"))
		((Real,radiusStdDev,0.3,,"Normalized standard deviation of generated sizes."))
		((Real,radiusMean,-1,,"Mean radius. If negative (default), autocomputed to as a function of box size and :yref:`TriaxialTest::numberOfGrains`"))
		((Real,sphereYoungModulus,15000000.0,,"Stiffness of spheres."))
		((Real,sphereKsDivKn,0.5,,"Ratio of shear vs. normal contact stiffness for spheres."))
		((Real,sphereFrictionDeg,18.0,,"Friction angle [°] of spheres assigned just before triaxial testing."))
		((Real,compactionFrictionDeg,sphereFrictionDeg,,"Friction angle [°] of spheres during compaction (different values result in different porosities)]. This value is overridden by :yref:`TriaxialTest::sphereFrictionDeg` before triaxial testing."))
		((Real,boxYoungModulus,15000000.0,,"Stiffness of boxes."))
		((Real,maxWallVelocity,10,,"max velocity of boundaries. Usually useless, but can help stabilizing the system in some cases."))
		((Real,boxKsDivKn,0.5,,"Ratio of shear vs. normal contact stiffness for boxes."))
		((Real,boxFrictionDeg,0.0,,"Friction angle [°] of boundaries contacts."))
		((Real,density,2600,,"density of spheres"))
		((Real,strainRate,0.1,,"Strain rate in triaxial loading."))
		((Real,defaultDt,-1,,"Max time-step. Used as initial value if defined. Latter adjusted by the time stepper."))
		((Real,dampingForce,0.2,,"Coefficient of Cundal-Non-Viscous damping (applied on on the 3 components of forces)"))
		((Real,dampingMomentum,0.2,,"Coefficient of Cundal-Non-Viscous damping (applied on on the 3 components of torques)"))
		((Real,StabilityCriterion,0.01,,"Value of unbalanced force for which the system is considered stable. Used in conditionals to switch between loading stages."))
		((Real,wallOversizeFactor,1.3,,"Make boundaries larger than the packing to make sure spheres don't go out during deformation."))
		((Real,sigmaIsoCompaction,50000,,"Confining stress during isotropic compaction."))
		((Real,sigmaLateralConfinement,50000,,"Lateral stress during triaxial loading. An isotropic unloading is performed if the value is not equal to :yref:`TriaxialTest::sigmaIsoCompaction`."))
		
		((int,timeStepUpdateInterval,50,,"interval for :yref:`GlobalStiffnessTimeStepper`"))
		((int,wallStiffnessUpdateInterval,10,,"interval for updating the stiffness of sample/boundaries contacts"))
		((int,radiusControlInterval,10,,"interval between size changes when growing spheres."))
		((int,numberOfGrains,400,,"Number of generated spheres."))
		((int,recordIntervalIter,20,,"interval between file outputs"))
		,
		/* init */
		,
		/* constructor */
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
		WallStressRecordFile = "./WallStresses"+Key;	
		gravity 		= Vector3r(0,-9.81,0);
		,
		//.def("setContactProperties",&TriaxialCompressionEngine::setContactProperties,"Assign a new friction angle (degrees) to dynamic bodies and relative interactions")
		 );
	DECLARE_LOGGER;
};

REGISTER_SERIALIZABLE(TriaxialTest);
