/*************************************************************************
*  Copyright (C) 2008 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include"StaticAttractionEngine.hpp"
#include<yade/core/Scene.hpp>


/*! \brief Compress stuff faster and better!
 *
 *   Faster compression is done by attracting everything to everything using
 *   classical gravitational equation, related to distance. And assuming that 
 *   everything has equal mass: 
 *                
 *                 D = 1/r^2,
 *
 *   where D is displacement delta, and r is radius. (Units intentionally do
 *   not add: [meter] = [1/m^2]. For unit correctness just multiply by [m^3] )
 *
 *   Forces are NOT used here. Only displacement value, which is normalized
 *   afterwards, which means that ALL values of D (for each body) are <= 1.
 *   Then all D (for all bodies) are multiplied by whatever returns the virtual
 *   function getMaxDisplacement(Scene*). Function doesItApplyToThisBody(Body*)
 *   is used to check if Body should be considered for displacement. This
 *   function can for example check if isDynamic.
 *
 *   This engine should be used at the end of the loop, because it directly
 *   changes position of bodies.
 *
 *   After the compression stage is over, engine should be disabled by 
 *   setting active=false;
 */
void StaticAttractionEngine::action()
{

	BodyContainer::iterator bi    = ncb->bodies->begin();
	BodyContainer::iterator biEnd = ncb->bodies->end();
	BodyContainer::iterator bi2;

	std::vector<Vector3r > moves;
	moves.clear();
	for( ; bi!=biEnd ; ++bi )
		if(doesItApplyToThisBody((*bi).get()))
		{
			Vector3r delta(0,0,0);
			Vector3r center_1 = (*bi)->physicalParameters->se3.position;
		
			bi2    = ncb->bodies->begin();
			for( ; bi2!=biEnd ; ++bi2 )
				if((*bi)->getId() != (*bi2)->getId() && doesItApplyToThisBody((*bi2).get()) ) // different bodies
				{
					Vector3r center_2 = (*bi2)->physicalParameters->se3.position;
					Vector3r dir=center_2-center_1;
					Real r = dir.Normalize(); // dir is unit vector, r is a distance
					delta += dir * 1/(r*r);
				}
			moves.push_back(delta);
		}

	Real maxl=0; // find maximum displacement
	for(unsigned int i = 0 ; i < moves.size() ; ++i )
		maxl = std::max(moves[i].Length(),maxl);

//	std::cerr << "maxl= " << maxl << "\n";
	for(unsigned int i = 0 ; i < moves.size() ; ++i ) // normalize, and multiply by planned maxium displacement
		moves[i] = getMaxDisplacement(ncb)*moves[i]/maxl;

	bi    = ncb->bodies->begin();
	for(unsigned int i = 0 ; bi!=biEnd ; ++bi )
		if(doesItApplyToThisBody((*bi).get()))
			(*bi)->physicalParameters->se3.position+=moves[i++];
}

YADE_PLUGIN((StaticAttractionEngine));

YADE_REQUIRE_FEATURE(PHYSPAR);

