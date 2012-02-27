/*************************************************************************
*  Copyright (C) 2009 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/core/InteractionSolver.hpp>

#include <set>
#include <boost/tuple/tuple.hpp>


class ElawSnowLayersDeformation : public InteractionSolver
{
/// Attributes
	public :
		int sdecGroupMask;
		Real creep_viscosity;
	
		ElawSnowLayersDeformation();
		void action();

	REGISTER_ATTRIBUTES(InteractionSolver,(sdecGroupMask)(creep_viscosity));
	REGISTER_CLASS_NAME(ElawSnowLayersDeformation);
	REGISTER_BASE_CLASS_NAME(InteractionSolver);
};

REGISTER_SERIALIZABLE(ElawSnowLayersDeformation);


