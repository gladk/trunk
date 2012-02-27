/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once


#include<yade/core/GeometricalModel.hpp>


class LineSegment : public GeometricalModel
{
	public :
		Real length;
		LineSegment ();
		virtual ~LineSegment ();
	
/// Serialization
	REGISTER_ATTRIBUTES(GeometricalModel,(length));
	REGISTER_CLASS_NAME(LineSegment);
	REGISTER_BASE_CLASS_NAME(GeometricalModel);
/// Indexable
	REGISTER_CLASS_INDEX(LineSegment,GeometricalModel);

};

REGISTER_SERIALIZABLE(LineSegment);


