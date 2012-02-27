/*************************************************************************
*  Copyright (C) 2006 by Luc Scholtes                                    *
*  luc.scholtes@hmg.inpg.fr                                              *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#ifndef SAMPLE_CAPILLARY_PRESSURE_ENGINE_HPP
#define SAMPLE_CAPILLARY_PRESSURE_ENGINE_HPP

#include<yade/core/PartialEngine.hpp>
#include<yade/lib-base/Math.hpp>
#include<yade/pkg-dem/TriaxialStressController.hpp>
#include <string>

/*! \brief Isotropic compression + suction variation test */

class Law2_ScGeom_CapillaryPhys_Capillarity;

class SampleCapillaryPressureEngine : public TriaxialStressController
{
	public :
		//! is isotropicInternalCompactionFinished?
		bool Phase1;
		std::string Phase1End;
		//! is this the beginning of the simulation, after reading the scene?
		bool firstRun;
		
		shared_ptr<Law2_ScGeom_CapillaryPhys_Capillarity>  capillaryCohesiveLaw;
		//Law2_ScGeom_CapillaryPhys_Capillarity* capillaryCohesiveLaw; // which one is right?
		
		virtual ~SampleCapillaryPressureEngine();
		void updateParameters();
		virtual void action();
		
	YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(SampleCapillaryPressureEngine,TriaxialStressController,"Rk: this engine has to be tested withthe new formalism. It produces the isotropic compaction of an assembly and allows one to controlled the capillary pressure inside (uses Law2_ScGeom_CapillaryPhys_Capillarity).",
		((Real,Pressure,0,,"Value of the capillary pressure Uc=Ugas-Uliquid (see Law2_ScGeom_CapillaryPhys_Capillarity). [Pa]"))
		((bool,pressureVariationActivated,1,,"Is the capillary pressure varying?"))
		((bool,fusionDetection,1,,"Is the detection of menisci overlapping activated?"))
		((bool,binaryFusion,1,,"If yes, capillary force are set to 0 when, at least, 1 overlap is detected for a meniscus. If no, capillary force is divided by the number of overlaps."))
		((Real,PressureVariation,0,,"Variation of the capillary pressure (each iteration). [Pa]"))
		((Real,UnbalancedForce,1,,"mean resultant forces divided by mean contact force"))
		((Real,StabilityCriterion,0.01,,"tolerance in terms of :yref:'TriaxialCompressionEngine::UnbalancedForce' to consider the packing as stable"))
		((Real,SigmaPrecision,0.001,,"tolerance in terms of mean stress to consider the packing as stable")),
		Phase1=false;
		Phase1End = "Compacted";
		firstRun=true;
		,)
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(SampleCapillaryPressureEngine);

#endif //  SAMPLECAPILLARYPRESSUREENGINE_HPP



