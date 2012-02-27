/*************************************************************************
*  Copyright (C) 2008 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/pkg-common/GLDrawFunctors.hpp>
#include<yade/pkg-common/Gl1_Sphere.hpp>


class Ef1_BssSnowGrain_glDraw : public GlShapeFunctor
{
//	Gl1_Sphere s;
	public :
		virtual void go(const shared_ptr<Shape>&, const shared_ptr<PhysicalParameters>&,bool,const GLViewInfo&);

	RENDERS(BssSnowGrain);
	REGISTER_CLASS_NAME(Ef1_BssSnowGrain_glDraw);
	REGISTER_BASE_CLASS_NAME(GlShapeFunctor);
};

REGISTER_SERIALIZABLE(Ef1_BssSnowGrain_glDraw);

