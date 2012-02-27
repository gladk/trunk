/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/core/GeometricalModel.hpp>
#include <Wm3Math.h>
#include<yade/lib-base/yadeWm3.hpp>

class Scene;

class Quadrilateral : public GeometricalModel
{
	public :
		int	 n1,n2,n3,n4;	// nodes of quadrilateral
		Vector3r i1,i2,i3,i4;	// initial coordinates FIXME - this Quadrilateral is used only 
					// for strain map display. some design changes are needed here anyway.

		Quadrilateral();
		Quadrilateral(int p1,int p2,int p3,int p4,Scene* mb);
		virtual ~Quadrilateral();

	REGISTER_ATTRIBUTES(GeometricalModel,(n1)(n2)(n3)(n4)(i1)(i2)(i3)(i4));
	REGISTER_CLASS_NAME(Quadrilateral);
	REGISTER_BASE_CLASS_NAME(GeometricalModel);
	REGISTER_CLASS_INDEX(Quadrilateral,GeometricalModel);
};

REGISTER_SERIALIZABLE(Quadrilateral);


