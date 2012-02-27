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


class ElasticCohesiveLaw : public InteractionSolver
{
	private :
		bool first; // FIXME - remove that!

	public :
		int sdecGroupMask;
		bool momentRotationLaw;
		
		ElasticCohesiveLaw();
		void action();

	REGISTER_ATTRIBUTES(InteractionSolver,(sdecGroupMask)(momentRotationLaw));

	REGISTER_CLASS_NAME(ElasticCohesiveLaw);
	REGISTER_BASE_CLASS_NAME(InteractionSolver);

};

REGISTER_SERIALIZABLE(ElasticCohesiveLaw);


