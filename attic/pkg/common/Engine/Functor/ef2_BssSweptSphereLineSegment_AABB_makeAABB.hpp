/*************************************************************************
*  Copyright (C) 2008 by Vincent Richefeu                                *
*  vincent.richefeu@hmg.inpg.fr                                          *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
 
#pragma once


#include<yade/pkg-common/BoundFunctor.hpp>


class ef2_BssSweptSphereLineSegment_AABB_makeAABB : public BoundFunctor
{
	public :
		void go(	  const shared_ptr<Shape>& cm
				, shared_ptr<Bound>& bv
				, const Se3r& se3
				, const Body*	);
                FUNCTOR2D(BssSweptSphereLineSegment,Aabb);
                REGISTER_CLASS_NAME(ef2_BssSweptSphereLineSegment_AABB_makeAABB);
	        REGISTER_BASE_CLASS_NAME(BoundFunctor);
};

REGISTER_SERIALIZABLE(ef2_BssSweptSphereLineSegment_AABB_makeAABB);


