/*************************************************************************
*  Copyright (C) 2008 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/core/PartialEngine.hpp>
#include<yade/pkg-common/StaticAttractionEngine.hpp>

class StaticSpheresAttractionEngine : public StaticAttractionEngine
{
	public:
		StaticSpheresAttractionEngine() : max_displacement(0) {};
		Real max_displacement;
	protected :
		virtual Real getMaxDisplacement(Scene*);
		virtual bool doesItApplyToThisBody(Body*);
	REGISTER_ATTRIBUTES(StaticAttractionEngine,/*(applies_to_what)(divide_radius)*/(max_displacement));
	REGISTER_CLASS_NAME(StaticSpheresAttractionEngine);
	REGISTER_BASE_CLASS_NAME(StaticAttractionEngine);
};

REGISTER_SERIALIZABLE(StaticSpheresAttractionEngine);


