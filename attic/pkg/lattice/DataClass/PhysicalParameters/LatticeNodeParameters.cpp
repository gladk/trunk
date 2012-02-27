/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "LatticeNodeParameters.hpp"

LatticeNodeParameters::LatticeNodeParameters() : PhysicalParameters()
{
	createIndex();
	
	countIncremental 	= 0;
	countStiffness 		= 0;
	displacementIncremental = Vector3r(0.0,0.0,0.0);
	displacementAlignmental 	= Vector3r(0.0,0.0,0.0);
}

LatticeNodeParameters::~LatticeNodeParameters()
{

}


YADE_PLUGIN((LatticeNodeParameters));

YADE_REQUIRE_FEATURE(PHYSPAR);

