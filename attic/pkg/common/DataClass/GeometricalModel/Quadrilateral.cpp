/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "Quadrilateral.hpp"
#include<yade/core/Scene.hpp>

YADE_REQUIRE_FEATURE(geometricalmodel);

Quadrilateral::Quadrilateral() : GeometricalModel()
{
	createIndex();
}

Quadrilateral::Quadrilateral(int p1,int p2,int p3,int p4,Scene* mb) : GeometricalModel()
{
	createIndex();
	n1 = p1;
	n2 = p2;
	n3 = p3;
	n4 = p4;
	#define Q(x) ( (*(mb->bodies))[x]->physicalParameters->se3.position )
	i1 = Q(n1);
	i2 = Q(n2);
	i3 = Q(n3);
	i4 = Q(n4);
	#undef Q
}

Quadrilateral::~Quadrilateral(){}

YADE_PLUGIN((Quadrilateral));


YADE_REQUIRE_FEATURE(PHYSPAR);

