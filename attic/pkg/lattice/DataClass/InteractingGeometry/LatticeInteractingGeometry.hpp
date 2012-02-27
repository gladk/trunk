/*************************************************************************
*  Copyright (C) 2007 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/pkg-common/SceneShape.hpp>

class LatticeInteractingGeometry : public SceneShape
{
	public :
		LatticeInteractingGeometry ();
		virtual ~LatticeInteractingGeometry ();

/// Serialization
	REGISTER_CLASS_NAME(LatticeInteractingGeometry);
	REGISTER_BASE_CLASS_NAME(SceneShape);
	
/// Indexable
	REGISTER_CLASS_INDEX(LatticeInteractingGeometry,SceneShape);
};

REGISTER_SERIALIZABLE(LatticeInteractingGeometry);


