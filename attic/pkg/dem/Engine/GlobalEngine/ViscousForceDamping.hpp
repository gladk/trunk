/*************************************************************************
*  Copyright (C) 2008 by Feng Chen  (fchen3@gmail.com)                   *
*  Department of Civil and Environmental Engineering                     *
*  223 Perkins Hall                                                      *
*  University of Tennessee, Knoxville, 37996                             *
*  http://fchen3.googlepages.com/home                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/core/InteractionSolver.hpp>

#include <set>
#include <boost/tuple/tuple.hpp>

class ViscousForceDamping : public InteractionSolver
{
/// Attributes
	public :
		
		Real betaNormal;
		Real betaShear;
		
		int sdecGroupMask;
		bool momentRotationLaw;
	
		ViscousForceDamping();
		void action(Body* body);

	REGISTER_ATTRIBUTES(InteractionSolver,(sdecGroupMask)(momentRotationLaw)(betaNormal)(betaShear));
	REGISTER_CLASS_NAME(ViscousForceDamping);
	REGISTER_BASE_CLASS_NAME(InteractionSolver);
};

REGISTER_SERIALIZABLE(ViscousForceDamping);


