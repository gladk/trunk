/*************************************************************************
*  Copyright (C) 2007 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/core/PartialEngine.hpp>
#include<Wm3Vector3.h>
#include<yade/lib-base/yadeWm3.hpp>

class DisplacementToForceEngine : public PartialEngine
{
	private :
		Real			target_length_sq,direction,old_direction,oscillations;
	public :
		Vector3r		targetForce,targetForceMask;
		Real			displacement;
		Vector3r		translationAxis;

		DisplacementToForceEngine();
		virtual ~DisplacementToForceEngine();

		virtual void action(Scene *);

	protected :
		virtual void postProcessAttributes(bool);
	REGISTER_ATTRIBUTES(PartialEngine,(displacement)(translationAxis)(targetForce)(targetForceMask));
	REGISTER_CLASS_NAME(DisplacementToForceEngine);
	REGISTER_BASE_CLASS_NAME(PartialEngine);
};

REGISTER_SERIALIZABLE(DisplacementToForceEngine);


