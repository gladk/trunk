/*************************************************************************
*  Copyright (C) 2007 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include"DisplacementToForceEngine.hpp"
#include<yade/core/Scene.hpp>
#include<yade/pkg-common/ParticleParameters.hpp>

DisplacementToForceEngine::DisplacementToForceEngine() : targetForce(Vector3r::Zero()), targetForceMask(Vector3r::Zero()) 
{
	direction=1.0;
	old_direction=1.0;
	target_length_sq=0;
	oscillations=1.0;
}

DisplacementToForceEngine::~DisplacementToForceEngine()
{
}

void DisplacementToForceEngine::postProcessAttributes(bool deserializing)
{
	if(deserializing)
	{
		translationAxis.Normalize();

		std::cerr << "displacement: " << displacement << "\n";
		std::cerr << "translationAxis: " << translationAxis<< "\n";
		std::cerr << "targetForce: " << targetForce<< "\n";
		std::cerr << "targetForceMask: " << targetForceMask<< "\n\n";
			
		target_length_sq = targetForce.SquaredLength();
	}
}



void DisplacementToForceEngine::action()
{
	shared_ptr<BodyContainer>& bodies = ncb->bodies;

	std::vector<int>::const_iterator ii = subscribedBodies.begin();
	std::vector<int>::const_iterator iiEnd = subscribedBodies.end();

	ncb->forces.sync();

	for(;ii!=iiEnd;++ii)
		if( bodies->exists(*ii) )
		{
			Vector3r current_force=ncb->forces.getForce(*ii);

			Real current_length_sq = 
				  (targetForceMask[0] != 0) ? current_force[0]*current_force[0]:0.0 
				+ (targetForceMask[1] != 0) ? current_force[1]*current_force[1]:0.0 
				+ (targetForceMask[2] != 0) ? current_force[2]*current_force[2]:0.0
				;
			old_direction=direction;
			direction = (current_force.Dot(targetForce) > 0) ? ( current_length_sq > target_length_sq ? -1.0 : 1.0 ) : 1.0 ;
			
			if(old_direction*direction < 0) oscillations+=1.0; else oscillations-=1.0;
			if(oscillations<1.0) oscillations = 1.0;
			//if(current_length_sq==0.0) direction=1.0;

			((*bodies)[*ii]->physicalParameters.get())->se3.position += displacement*translationAxis*direction/oscillations;
		}
		
}

YADE_PLUGIN((DisplacementToForceEngine));

YADE_REQUIRE_FEATURE(PHYSPAR);

