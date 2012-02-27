/*************************************************************************
*  Copyright (C) 2008 by Vincent Richefeu                                *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/


#include "ef2_BshTube_BssSweptSphereLineSegment_makeBssSweptSphereLineSegment.hpp"
#include<yade/pkg-common/BssSweptSphereLineSegment.hpp>
#include<yade/pkg-common/BshTube.hpp>
#include<yade/pkg-common/BcpConnection.hpp>
YADE_REQUIRE_FEATURE(geometricalmodel);

void ef2_BshTube_BssSweptSphereLineSegment_makeBssSweptSphereLineSegment::go(	const shared_ptr<GeometricalModel>& gm,
				shared_ptr<Shape>& ig,
				const Se3r& se3,
				const Body* body)
{
  
  BshTube* tube = static_cast<BshTube*>(gm.get());
  if(ig == 0)
    ig = boost::shared_ptr<Shape>(new BssSweptSphereLineSegment);

  BssSweptSphereLineSegment* SSLS = static_cast<BssSweptSphereLineSegment*>(ig.get());
  
  SSLS->radius = tube->radius;
  SSLS->length = 2.0 * tube->half_height;
  
  BcpConnection* bc = static_cast<BcpConnection*>(body->physicalParameters.get());

  SSLS->position = (*(Omega::instance().getScene()->bodies))[bc->id1]->physicalParameters->se3.position;
  Vector3r len = (*(Omega::instance().getScene()->bodies))[bc->id2]->physicalParameters->se3.position - SSLS->position;
  len.Normalize();
  SSLS->orientation = len;
  
}
	
YADE_PLUGIN((ef2_BshTube_BssSweptSphereLineSegment_makeBssSweptSphereLineSegment));

YADE_REQUIRE_FEATURE(PHYSPAR);

