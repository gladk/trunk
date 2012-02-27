/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/pkg-common/GLDrawFunctors.hpp>

class GLDrawLatticeBeamState : public GlStateFunctor
{
	public :
		virtual void go(const shared_ptr<PhysicalParameters>&);

	RENDERS(LatticeBeamParameters);
	REGISTER_CLASS_NAME(GLDrawLatticeBeamState);
	REGISTER_BASE_CLASS_NAME(GlStateFunctor);

};

REGISTER_SERIALIZABLE(GLDrawLatticeBeamState);


