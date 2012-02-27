/*************************************************************************
*  Copyright (C) 2007 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "LatticeInteractingGeometry.hpp"


LatticeInteractingGeometry::LatticeInteractingGeometry () : SceneShape()
{
	createIndex();
}

LatticeInteractingGeometry::~LatticeInteractingGeometry ()
{
}


YADE_PLUGIN((LatticeInteractingGeometry));

YADE_REQUIRE_FEATURE(PHYSPAR);

