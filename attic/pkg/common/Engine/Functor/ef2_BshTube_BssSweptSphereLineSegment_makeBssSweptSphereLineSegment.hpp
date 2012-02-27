/*************************************************************************
*  Copyright (C) 2008 by Vincent Richefeu                                *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/pkg-common/InteractingGeometryEngineUnit.hpp>

class ef2_BshTube_BssSweptSphereLineSegment_makeBssSweptSphereLineSegment : public InteractingGeometryEngineUnit
{
	public :
		void go(	  const shared_ptr<GeometricalModel>& gm
				, shared_ptr<Shape>& ig
				, const Se3r& se3
				, const Body*	);

        FUNCTOR2D(BshTube,BssSweptSphereLineSegment);
        REGISTER_CLASS_NAME(ef2_BshTube_BssSweptSphereLineSegment_makeBssSweptSphereLineSegment);
	REGISTER_BASE_CLASS_NAME(InteractingGeometryEngineUnit);
	
        DEFINE_FUNCTOR_ORDER_2D(BshTube,BssSweptSphereLineSegment);
};

REGISTER_SERIALIZABLE(ef2_BshTube_BssSweptSphereLineSegment_makeBssSweptSphereLineSegment);



