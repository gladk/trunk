/*************************************************************************
*  Copyright (C) 2008 by Vincent Richefeu                                *
*  vincent.richefeu@hmg.inpg.fr                                          *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "BssSweptSphereLineSegment.hpp"


BssSweptSphereLineSegment::BssSweptSphereLineSegment () : Shape()
{		
	createIndex();
        position = orientation = Vector3r(0,0,0);
        radius = length = 0.0;
}

BssSweptSphereLineSegment::~BssSweptSphereLineSegment ()
{		
}

YADE_PLUGIN((BssSweptSphereLineSegment));

YADE_REQUIRE_FEATURE(PHYSPAR);

