/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "LatticeBeamAngularSpring.hpp"


LatticeBeamAngularSpring::LatticeBeamAngularSpring() : InteractionPhysics()
{

	initialPlaneAngle	= 0;
//	planeAngle		= 0;
	lastCrossProduct	= Vector3r(0,0,0); // stupid value, zero is useless...
	initialOffPlaneAngle1	= 0;
	initialOffPlaneAngle2	= 0;
	lastOffPlaneAngleDifference1	= 0;
	lastOffPlaneAngleDifference2	= 0;
	swirl1			= 0;
	swirl2			= 0;
//	offPlaneAngle		= 0;
	planeSwap180		= false;

	createIndex();
}


LatticeBeamAngularSpring::~LatticeBeamAngularSpring()
{

}


YADE_PLUGIN((LatticeBeamAngularSpring));

YADE_REQUIRE_FEATURE(PHYSPAR);

