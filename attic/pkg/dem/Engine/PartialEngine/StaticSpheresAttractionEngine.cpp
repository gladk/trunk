/*************************************************************************
*  Copyright (C) 2008 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include <yade/pkg-dem/StaticSpheresAttractionEngine.hpp>
#include <yade/pkg-dem/FrictPhys.hpp>
#include <yade/pkg-common/Sphere.hpp>
#include <yade/core/Scene.hpp>

Real StaticSpheresAttractionEngine::getMaxDisplacement()
{
	// calculation based on max normal force - FIXME not working satifactionary, currently
	/*
	Real penetration=0;
	InteractionContainer::iterator in     = ncb->interactions->begin();
	InteractionContainer::iterator in_end = ncb->interactions->end();
	for( ; in!=in_end ; ++in )
		if(	   doesItApplyToThisBody((*(ncb->bodies))[ (*in)->getId1() ].get()) 
			&& doesItApplyToThisBody((*(ncb->bodies))[ (*in)->getId2() ].get()))
		{ // FIXME - it should work with more than FrictPhys
			FrictPhys* el = dynamic_cast<FrictPhys*>((*in)->interactionPhysics.get());
			if(el)
				penetration = std::max(el->normalForce.Length() / el->initialKn , penetration);
		}

	max_displacement = penetration / 10.0;
	*/

	if(max_displacement == 0) // OK, so try radius based approach
	{
		// simple calculation based on max radius
		Real max_radius=0;
		BodyContainer::iterator bi    = ncb->bodies->begin();
		BodyContainer::iterator biEnd = ncb->bodies->end();
		for( ; bi!=biEnd ; ++bi )
			if(doesItApplyToThisBody((*bi).get()))
				max_radius = std::max(static_cast<Sphere*>((*bi)->shape.get())->radius , max_radius);
		max_displacement = max_radius / 50;
	}
	return max_displacement;
};

bool StaticSpheresAttractionEngine::doesItApplyToThisBody(Body* b)
{
	if(b->isDynamic && b->shape->getClassName()=="Sphere")
		return true; 
	return false;
};


YADE_PLUGIN((StaticSpheresAttractionEngine));

YADE_REQUIRE_FEATURE(PHYSPAR);

