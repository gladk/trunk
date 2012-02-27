/*************************************************************************
*  Copyright (C) 2008 by Vincent Richefeu                                *
*  vincent.richefeu@hmg.inpg.fr                                          *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/pkg-common/GLDrawFunctors.hpp>

class GLDrawBssSweptSphereLineSegment : public GlShapeFunctor
{	
	private :
		static bool first;
		static int  glWiredList;
		static int  glList;
	
	public :
                GLDrawBssSweptSphereLineSegment();
		virtual void go(const shared_ptr<Shape>&, const shared_ptr<PhysicalParameters>&,bool,const GLViewInfo&);
/// Serialization
	protected :
		virtual void postProcessAttributes(bool deserializing){if(deserializing){first=true;};};

        RENDERS(BssSweptSphereLineSegment);
        REGISTER_CLASS_NAME(GLDrawBssSweptSphereLineSegment);
	REGISTER_BASE_CLASS_NAME(GlShapeFunctor);
};

REGISTER_SERIALIZABLE(GLDrawBssSweptSphereLineSegment);


