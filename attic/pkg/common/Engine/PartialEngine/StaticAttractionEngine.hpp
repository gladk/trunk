/*************************************************************************
*  Copyright (C) 2008 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/core/PartialEngine.hpp>
#include <Wm3Vector3.h>
#include<yade/lib-base/yadeWm3.hpp>

class StaticAttractionEngine : public PartialEngine
{
	public :
		void action(Scene *);
		bool active;
		StaticAttractionEngine() : active(true) {};
		virtual bool isActivated(Scene*) { return active; };

	protected :
		virtual Real getMaxDisplacement(Scene*) {throw "do not use StaticAttractionEngine directly, use derived classes instead.";};
		virtual bool doesItApplyToThisBody(Body*) {throw "do not use StaticAttractionEngine directly, use derived classes instead.";};
	REGISTER_ATTRIBUTES(PartialEngine,(active));
	REGISTER_CLASS_NAME(StaticAttractionEngine);
	REGISTER_BASE_CLASS_NAME(PartialEngine);
};

REGISTER_SERIALIZABLE(StaticAttractionEngine);


