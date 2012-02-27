/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/core/PhysicalParameters.hpp>
#include <vector>

class LatticeBeamParameters : public PhysicalParameters 
{
	private:
		long int                 lastIter_;
		Real                     strain_;
	public :
		// element
		int                      id1
					,id2;
		// state		
		Real  			 initialLength
					,length;

		Vector3r                 direction		// is a unit vector
					,otherDirection;	// is a unit vector too
		Vector3r		 shearing_strain;

		Se3r                     se3Displacement;

		// parameters
		Real                     criticalTensileStrain
					,criticalCompressiveStrain
	
					,longitudalStiffness
					,bendingStiffness
					,torsionalStiffness;


		// where does it belong, what is it?
		Real 			 count
					,torsionAngle;
		Quaternionr		 bendingRotation
					,torsionalRotation;

		Real strain();
		// must go to derived class
		//Real                     nonLocalStrain
		//                        ,nonLocalDivisor;

		LatticeBeamParameters();
		virtual ~LatticeBeamParameters();

/// Serialization
	protected :
		virtual void postProcessAttributes(bool);
	REGISTER_ATTRIBUTES(PhysicalParameters,
		(id1)
		(id2)
		
		(initialLength)
		(length) // FIXME - can be calculated after deserialization
		
		(direction) // FIXME - can be calculated after deserialization
		(otherDirection) // this must be saved!
		
		(criticalTensileStrain)
		(criticalCompressiveStrain)

		(longitudalStiffness)
		(bendingStiffness)
		(torsionalStiffness)

		(se3Displacement)
	);
	REGISTER_CLASS_NAME(LatticeBeamParameters);
	REGISTER_BASE_CLASS_NAME(PhysicalParameters);
	
/// Indexable
	REGISTER_CLASS_INDEX(LatticeBeamParameters,PhysicalParameters);

};

REGISTER_SERIALIZABLE(LatticeBeamParameters);


