/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/core/InteractionPhysics.hpp>
#include <vector>

class NonLocalDependency : public InteractionPhysics
{
	public :
		Real			gaussValue;
		
		NonLocalDependency();
		virtual ~NonLocalDependency();

/// Serialization
	REGISTER_ATTRIBUTES(InteractionPhysics,(gaussValue));
	REGISTER_CLASS_NAME(NonLocalDependency);
	REGISTER_BASE_CLASS_NAME(InteractionPhysics);
/// Indexable
	REGISTER_CLASS_INDEX(NonLocalDependency,InteractionPhysics);
};
REGISTER_SERIALIZABLE(NonLocalDependency);


