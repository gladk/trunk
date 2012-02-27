/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/core/InteractionPhysics.hpp>
#include <vector>

class LatticeBeamAngularSpring : public InteractionPhysics
{
	public :
		Vector3r		 lastCrossProduct;
		bool			 planeSwap180;
		Real			 initialPlaneAngle  // 2D only for now
			 		//,planeAngle
					,initialOffPlaneAngle1
					,initialOffPlaneAngle2
					,lastOffPlaneAngleDifference1
					,lastOffPlaneAngleDifference2
					//,offPlaneAngle
					;
		int			 swirl1
					,swirl2;
		
		LatticeBeamAngularSpring();
		virtual ~LatticeBeamAngularSpring();

/// Serialization
	REGISTER_ATTRIBUTES(InteractionPhysics,
		(initialPlaneAngle)
		// (planeAngle) // FIXME - can be calculated after deserialization
		(lastCrossProduct)
		(planeSwap180)
		(initialOffPlaneAngle1)
		(initialOffPlaneAngle2)
		(lastOffPlaneAngleDifference1)
		(lastOffPlaneAngleDifference2)
		(swirl1)
		(swirl2)
		//	(offPlaneAngle) // FIXME - can be calculated after deserialization
	);
	REGISTER_CLASS_NAME(LatticeBeamAngularSpring);
	REGISTER_BASE_CLASS_NAME(InteractionPhysics);
	
/// Indexable
	REGISTER_CLASS_INDEX(LatticeBeamAngularSpring,InteractionPhysics);

};

REGISTER_SERIALIZABLE(LatticeBeamAngularSpring);


