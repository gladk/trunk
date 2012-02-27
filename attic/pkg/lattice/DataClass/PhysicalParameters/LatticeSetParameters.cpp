/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "LatticeSetParameters.hpp"


LatticeSetParameters::LatticeSetParameters() : PhysicalParameters()
{
        createIndex();
        nodeGroupMask = 1;
        beamGroupMask = 2;
        total = 0;
        nonl = 0;
        range=0;
        useBendTensileSoftening = false;
        useStiffnessSoftening   = false;
}


LatticeSetParameters::~LatticeSetParameters()
{

}


 
YADE_PLUGIN((LatticeSetParameters));

YADE_REQUIRE_FEATURE(PHYSPAR);

