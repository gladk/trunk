/*************************************************************************
*  Copyright (C) 2004 by Olivier Galizzi                                 *
*  olivier.galizzi@imag.fr                                               *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/core/InteractionSolver.hpp>

#include <set>
#include <boost/tuple/tuple.hpp>

class VolumicContactLaw : public InteractionSolver
{
/// Attributes
	public :
		int sdecGroupMask;
		bool momentRotationLaw;
		vector<double> compacite_init;
		Real rigidite_max_enreg;
		Real max_force_enreg;
		Real compacity_init;
		Real compacity;
		Real compacity_max;
		Real interpenetration_max;
		Real delta_interpenetration_max;


	
		VolumicContactLaw();
		void action();
		void speedTest();
	REGISTER_ATTRIBUTES(InteractionSolver,(sdecGroupMask)(momentRotationLaw)/*(compacite_init)*/);
	REGISTER_CLASS_NAME(VolumicContactLaw);
	REGISTER_BASE_CLASS_NAME(InteractionSolver);
};

REGISTER_SERIALIZABLE(VolumicContactLaw);


