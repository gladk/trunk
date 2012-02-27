/*************************************************************************
*  Copyright (C) 2006 by Bruno Chareyre		                         *
*  bruno.chareyre@hmg.inpg.fr                                            *
*  Copyright (C) 2006 by luc Scholtes                                    *
*  luc.scholtes@hmg.inpg.fr                                              *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/core/FileGenerator.hpp>
#include<yade/lib-base/Math.hpp>

class ForceRecorder;
class AveragePositionRecorder;
class VelocityRecorder;
class TriaxialStressController;
class TriaxialCompressionEngine;
class TriaxialStateRecorder;
class CapillaryStressRecorder;
class GlobalStiffnessTimeStepper;

/*! \brief Triaxial test on unsaturated sphere packings

	This preprocessor is a variant of TriaxialTest, including the model of capillary forces developed as part of the PhD of Luc Scholtès. See the documentation of Law2_ScGeom_CapillaryPhys_Capillarity or the main page https://yade-dem.org/wiki/CapillaryTriaxialTest, for more details.
	Results obtained with this preprocessor were reported for instance in "Scholtes et al. Micromechanics of granular materials with capillary effects. International Journal of Engineering Science 2009,(47)1, 64-75."
	
 */

class CapillaryTriaxialTest : public FileGenerator
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
		shared_ptr<CapillaryStressRecorder> capillaryStressRecorder;
		shared_ptr<GlobalStiffnessTimeStepper> globalStiffnessTimeStepper;
			
		void createBox(shared_ptr<Body>& body, Vector3r position, Vector3r extents,bool wire);
		void createSphere(shared_ptr<Body>& body, Vector3r position, Real radius,bool big,bool dynamic);
		void createActors(shared_ptr<Scene>& scene);
		void positionRootBody(shared_ptr<Scene>& scene);
		typedef pair<Vector3r, Real> BasicSphere;	
	public : 
		~CapillaryTriaxialTest ();
		bool generate(std::string& message);
		
		YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(
		CapillaryTriaxialTest,FileGenerator,"This preprocessor is a variant of TriaxialTest, including the model of capillary forces developed as part of the PhD of Luc Scholtès. See the documentation of Law2_ScGeom_CapillaryPhys_Capillarity or the main page https://yade-dem.org/wiki/CapillaryTriaxialTest, for more details.\n\n Results obtained with this preprocessor were reported for instance in 'Scholtes et al. Micromechanics of granular materials with capillary effects. International Journal of Engineering Science 2009,(47)1, 64-75.'"
		,
   		((Vector3r,lowerCorner,Vector3r(0,0,0),,"Lower corner of the box."))
		((Vector3r,upperCorner,Vector3r(1,1,1),,"Upper corner of the box."))
		((string,importFilename,"",,"File with positions and sizes of spheres."))
		((string,Key,"",,"A code that is added to output filenames."))
		((string,fixedBoxDims,"",,"string that contains some subset (max. 2) of {'x','y','z'} ; contains axes will have box dimension hardcoded, even if box is scaled as mean_radius is prescribed: scaling will be applied on the rest."))
		((Real,CapillaryPressure,0,,"Define succion in the packing [Pa]. This is the value used in the capillary model."))
		((bool,water,true,,"activate capillary model"))
		((bool,fusionDetection,false,,"test overlaps between liquid bridges on modify forces if overlaps exist"))
		((bool,binaryFusion,true,,"Defines how overlapping bridges affect the capillary forces (see :yref:`CapillaryTriaxialTest::fusionDetection`). If binary=true, the force is null as soon as there is an overlap detected, if not, the force is divided by the number of overlaps."))
		((string,WallStressRecordFile,"./WallStressesWater"+Key,,""))
		((string,capillaryStressRecordFile,"./capStresses"+Key,,""))
		((string,contactStressRecordFile,"./contStresses"+Key,,""))
		((bool,internalCompaction,false,,"flag for choosing between moving boundaries or increasing particles sizes during the compaction stage."))
		((bool,biaxial2dTest,false,,"FIXME : what is that?"))
		((bool,fixedPoroCompaction,false,,"flag to choose an isotropic compaction until a fixed porosity choosing a same translation speed for the six walls"))
		((bool,autoCompressionActivation,true,,"Do we just want to generate a stable packing under isotropic pressure (false) or do we want the triaxial loading to start automatically right after compaction stage (true)?"))
		((bool,autoUnload,true,,"auto adjust the isotropic stress state from :yref:`TriaxialTest::sigmaIsoCompaction` to :yref:`TriaxialTest::sigmaLateralConfinement` if they have different values. See docs for :yref:`TriaxialCompressionEngine::autoUnload`"))
		((bool,autoStopSimulation,false,,"freeze the simulation when conditions are reached (don't activate this if you want to be able to run/stop from Qt GUI)"))
		((bool,noFiles,false,,"Do not create any files during run (.xml, .spheres, wall stress records)"))
		((bool,facetWalls,false,,"Use facets for boundaries (not tested)"))
		((bool,wallWalls,false,,"Use walls for boundaries (not tested)"))
		((bool,boxWalls,true,,"Use boxes for boundaries (recommended)."))		
		((Real,fixedPorosity,1,,"FIXME : what is that?"))
		((Real,thickness,0.001,,"thickness of boundaries. It is arbitrary and should have no effect"))
		((Real,maxMultiplier,1.01,,"max multiplier of diameters during internal compaction (initial fast increase)"))
		((Real,finalMaxMultiplier,1.001,,"max multiplier of diameters during internal compaction (secondary precise adjustment)"))
		((Real,Rdispersion,0.3,,"Normalized standard deviation of generated sizes."))
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
		((Real,strainRate,1,,"Strain rate in triaxial loading."))
		((Real,defaultDt,0.0001,,"Max time-step. Used as initial value if defined. Latter adjusted by the time stepper."))
		((Real,dampingForce,0.2,,"Coefficient of Cundal-Non-Viscous damping (applied on on the 3 components of forces)"))
		((Real,dampingMomentum,0.2,,"Coefficient of Cundal-Non-Viscous damping (applied on on the 3 components of torques)"))
		((Real,StabilityCriterion,0.01,,"Value of unbalanced force for which the system is considered stable. Used in conditionals to switch between loading stages."))
		((Real,wallOversizeFactor,1.3,,"Make boundaries larger than the packing to make sure spheres don't go out during deformation."))
		((Real,sigmaIsoCompaction,50000,,"Confining stress during isotropic compaction."))
		((Real,sigmaLateralConfinement,50000,,"Lateral stress during triaxial loading. An isotropic unloading is performed if the value is not equal to :yref:`CapillaryTriaxialTest::SigmaIsoCompaction`."))
		
		((int,timeStepUpdateInterval,50,,"interval for :yref:`GlobalStiffnessTimeStepper`"))
		((int,timeStepOutputInterval,50,,"interval for outputing general information on the simulation (stress,unbalanced force,...)"))
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
// 	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(CapillaryTriaxialTest);


