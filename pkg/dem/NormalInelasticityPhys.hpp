/*************************************************************************
*  Copyright (C) 2008 by Jérôme DURIEZ                                   *
*  duriez@geo.hmg.inpg.fr                                                *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/pkg-dem/FrictPhys.hpp>


class NormalInelasticityPhys : public FrictPhys
{
	public :
		virtual ~NormalInelasticityPhys();

	YADE_CLASS_BASE_DOC_ATTRS_CTOR(NormalInelasticityPhys,FrictPhys,
				"Physics (of interaction) for using :yref:`Law2_ScGeom6D_NormalInelasticityPhys_NormalInelasticity` : with inelastic unloadings",
				((Real,unMax,0.0,,"the maximum value of penetration depth of the history of this interaction"))
				((Real,previousun,0.0,,"the value of this un at the last time step"))
				((Real,previousFn,0.0,,"the value of the normal force at the last time step"))
				((Real,forMaxMoment,1.0,,"parameter stored for each interaction, and allowing to compute the maximum value of the exchanged torque : TorqueMax= forMaxMoment * NormalForce"))
				((Real,kr,0.0,,"the rolling stiffness of the interaction"))
				((Real,knLower,0.0,,"the stifness corresponding to a virgin load for example"))
				// internal attributes
				((Vector3r,moment_twist,Vector3r(0,0,0),(Attr::noSave | Attr::readonly),"Twist moment. Defined here, being initialized as it should be, to be used in :yref:`Law2_ScGeom6D_NormalInelasticityPhys_NormalInelasticity`"))
				((Vector3r,moment_bending,Vector3r(0,0,0),(Attr::noSave | Attr::readonly),"Bending moment. Defined here, being initialized as it should be, to be used in :yref:`Law2_ScGeom6D_NormalInelasticityPhys_NormalInelasticity`"))
				,
				createIndex();
				);
	REGISTER_CLASS_INDEX(NormalInelasticityPhys,FrictPhys);
};

REGISTER_SERIALIZABLE(NormalInelasticityPhys);


