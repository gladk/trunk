/*************************************************************************
*  Copyright (C) 2008 by Vincent Richefeu                                *
*  vincent.richefeu@hmg.inpg.fr                                          *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/


#pragma once


#include<yade/core/GeometricalModel.hpp>


class BshTube : public GeometricalModel
{
	public :
		Real  radius;
                Real  half_height;

		BshTube ();
		virtual ~BshTube ();
	
	REGISTER_CLASS_NAME(BshTube);
	REGISTER_BASE_CLASS_NAME(GeometricalModel);
	REGISTER_ATTRIBUTES(GeometricalModel,(radius)(half_height));
	REGISTER_CLASS_INDEX(BshTube,GeometricalModel);
};

REGISTER_SERIALIZABLE(BshTube);


