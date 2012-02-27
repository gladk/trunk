/*************************************************************************
*  Copyright (C) 2009 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
 
#pragma once

#include<yade/pkg-common/BoundFunctor.hpp>

class Ef2_BssSnowGrain_AABB_makeAABB : public BoundFunctor
{
	public :
		Ef2_BssSnowGrain_AABB_makeAABB() {};
		void go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se3r& se3, const Body*);

	FUNCTOR2D(BssSnowGrain,Aabb);
	REGISTER_CLASS_NAME(Ef2_BssSnowGrain_AABB_makeAABB);
	REGISTER_BASE_CLASS_NAME(BoundFunctor);
};

REGISTER_SERIALIZABLE(Ef2_BssSnowGrain_AABB_makeAABB);

