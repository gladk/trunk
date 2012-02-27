/*************************************************************************
*  Copyright (C) 2008 by Vincent Richefeu                                *
*  vincent.richefeu@hmg.inpg.fr                                          *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/


#include "BshTube.hpp"

YADE_REQUIRE_FEATURE(geometricalmodel);

BshTube::BshTube () : GeometricalModel()
{		
	createIndex();
}


BshTube::~BshTube ()
{		
}

YADE_PLUGIN((BshTube));


YADE_REQUIRE_FEATURE(PHYSPAR);

