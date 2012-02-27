/*************************************************************************
*  Copyright (C) 2008 by Sergei Dorofeenko				 *
*  sega@users.berlios.de                                                 *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
 
#include <yade/pkg-common/Facet.hpp>
#include "Bo1_Facet_Aabb.hpp"
#include <yade/pkg-common/Aabb.hpp>

void Bo1_Facet_Aabb::go(	  const shared_ptr<Shape>& cm
				, shared_ptr<Bound>& bv
				, const Se3r& se3
				, const Body*	)
{
	Aabb* aabb = static_cast<Aabb*>(bv.get());
	Facet* facet = static_cast<Facet*>(cm.get());
	const Vector3r& O = se3.position;
	Matrix3r facetAxisT=se3.orientation.toRotationMatrix();
	const vector<Vector3r>& vertices=facet->vertices;
	if(!scene->isPeriodic){
		aabb->min=aabb->max = O + facetAxisT * vertices[0];
		for (int i=1;i<3;++i)
		{
			Vector3r v = O + facetAxisT * vertices[i];
			aabb->min = componentMinVector( aabb->min, v);
			aabb->max = componentMaxVector( aabb->max, v);
		}
	} else {
		Real inf=std::numeric_limits<Real>::infinity();
		aabb->min=Vector3r(inf,inf,inf); aabb->max=Vector3r(-inf,-inf,-inf);
		for(int i=0; i<3; i++){
			Vector3r v=scene->cell->unshearPt(O+facetAxisT*vertices[i]);
			aabb->min=componentMinVector(aabb->min,v);
			aabb->max=componentMaxVector(aabb->max,v);
		}
	}
}
	
YADE_PLUGIN((Bo1_Facet_Aabb));
