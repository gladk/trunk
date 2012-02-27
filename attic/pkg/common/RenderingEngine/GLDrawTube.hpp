/*************************************************************************
*  Copyright (C) 2008 by Vincent Richefeu                                *
*  vincent.richefeu@hmg.inpg.fr                                          *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/pkg-common/GLDrawFunctors.hpp>


class GLDrawTube : public GlGeometricalModelFunctor
{
	private :
		static bool first;
		static int  glWiredTubeList;
		static int  glTubeList;
                
	public :
		virtual void go(const shared_ptr<GeometricalModel>&, const shared_ptr<PhysicalParameters>&,bool);

                GLDrawTube() { first = true; }
                
	RENDERS(BshTube);
	REGISTER_CLASS_NAME(GLDrawTube);
	REGISTER_BASE_CLASS_NAME(GlGeometricalModelFunctor);
};

REGISTER_SERIALIZABLE(GLDrawTube);


