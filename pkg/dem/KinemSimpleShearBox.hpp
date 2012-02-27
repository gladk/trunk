/*************************************************************************
*  Copyright (C) 2010 by Jerome Duriez                                   *
*  jerome.duriez@hmg.inpg.fr                                             *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/pkg-common/BoundaryController.hpp>
#include<yade/core/Body.hpp>
#include<yade/lib-base/Math.hpp>
#include <yade/core/Scene.hpp>




class KinemSimpleShearBox : public BoundaryController
{

	protected :
		Real Scontact // computed by computeScontact()
		     ,dt // value of the time step, set in setBoxes_Dt
		     ,stiffness	// the normal stifness on the sample acting below the upper plate. Computed by computeStiffness
		     ,dalpha	// the increment over alpha, due to vertical displacement of upper box
		      ,deltaH	// the vertical increment of displacement to allow on the upper plate on this time step to verify either the constant normal stress or the constant normal stifness. Computed here by computeDY(..)
			;
			
		shared_ptr<Body> leftbox;
		shared_ptr<Body> rightbox;
		shared_ptr<Body> frontbox;
		shared_ptr<Body> backbox;
		shared_ptr<Body> topbox;
		shared_ptr<Body> boxbas;

		void 	computeAlpha()
			,computeScontact()
			,stopMovement() // to cancel all velocities when end of the loading is reached
			,letMove(Real dgamma,Real dH) //dgamma and dH being computed by different ways in the different Kinem... Engines
			,computeStiffness() // computes the stiffness of the contact sample - upper side. Useful for CNL and CNS loads
			,computeDY(Real KnC) // computes vertical displacement to perform to obey to the stiffness of the shear wanted : KnC. Useful for CNL (for which KnC = 0) and CNS loads
			,setBoxes_Dt()
			;
	
	YADE_CLASS_BASE_DOC_ATTRS_CTOR(KinemSimpleShearBox,BoundaryController,
			 "This class is supposed to be a mother class for all Engines performing loadings on the simple shear box of :yref:`SimpleShear`. It is not intended to be used by itself, but its declaration and implentation will thus contain all what is useful for all these Engines. The script simpleShear.py illustrates the use of the various corresponding Engines.",
			((Real,alpha,Mathr::PI/2.0,,"the angle from the lower box to the left box (trigo wise). Measured by this Engine, not to be changed by the user."))
			((std::vector<Real>,temoin_save,,,"vector (same length as 'gamma_save' for ex), with 0 or 1 depending whether the save for the corresponding value of gamma has been done (1) or not (0). Not to be changed by the user."))
			((Body::id_t,id_topbox,3,,"the id of the upper wall"))
			((Body::id_t,id_boxbas,1,,"the id of the lower wall"))
			((Body::id_t,id_boxleft,0,,"the id of the left wall"))
			((Body::id_t,id_boxright,2,,"the id of the right wall"))
			((Body::id_t,id_boxfront,5,,"the id of the wall in front of the sample"))
			((Body::id_t,id_boxback,4,,"the id of the wall at the back of the sample"))
			((Real,max_vel,1.0,,"to limit the speed of the vertical displacements done to control $\\sigma$ (CNL or CNS cases) [$m/s$]"))
			((Real,wallDamping,0.2,,"the vertical displacements done to to control $\\sigma$ (CNL or CNS cases) are in fact damped, through this wallDamping"))
			((bool,firstRun,true,,"boolean set to false as soon as the engine has done its job one time : useful to know if initial height of, and normal force sustained by, the upper box are known or not (and thus if they have to be initialized). Not to be changed by the user."))
			((Real,f0,0.0,,"the (vertical) force acting on the upper plate on the very first time step (determined by the Engine). Controls of the loadings in case of :yref:`KinemCNSEngine` or :yref:`KinemCNLEngine` will be done according to this initial value [$N$]. Not to be changed by the user.]"))
			((Real,y0,0.0,,"the height of the upper plate at the very first time step : the engine finds its value [$m$]. Not to be changed by the user."))
			((bool,LOG,false,,"boolean controling the output of messages on the screen")) //FIXME : surely something better to use here
			((string,Key,"",,"string to add at the names of the saved files")),
/*			leftbox = Body::byId(id_boxleft,scene);
			rightbox = Body::byId(id_boxright,scene);
			frontbox = Body::byId(id_boxfront,scene);
			backbox = Body::byId(id_boxback,scene);
			topbox = Body::byId(id_topbox,scene);
			boxbas = Body::byId(id_boxbas,scene);*/
			Scontact = 0.0;
			dt = 0.0;
			stiffness=0.0;
			dalpha=0.0;
			deltaH=0.0;
					 );
};

REGISTER_SERIALIZABLE(KinemSimpleShearBox);


