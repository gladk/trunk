/*************************************************************************
*  Copyright (C) 2008 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
 
#include"Ef2_BssSnowGrain_AABB_makeAABB.hpp"
#include<yade/pkg-common/Sphere.hpp>
#include<yade/pkg-snow/BssSnowGrain.hpp>
#include<yade/pkg-common/Aabb.hpp>
YADE_REQUIRE_FEATURE(geometricalmodel);
void Ef2_BssSnowGrain_AABB_makeAABB::go(const shared_ptr<Shape>& cm, shared_ptr<Bound>& bv, const Se3r& se3, const Body* b)
{
	BssSnowGrain* bss = static_cast<BssSnowGrain*>(cm.get());
	Aabb* aabb = static_cast<Aabb*>(bv.get());

/*
	public:
		Vector3r center,c_axis;
		Vector3r start,end;
		Vector3r color;
		int selection;
		std::vector<std::vector<Vector3r> > slices;
		Real layer_distance;

	m_copy.center   =grain->center;
	m_copy.c_axis   =grain->c_axis;
	m_copy.start    =grain->start;
	m_copy.end      =grain->end;
	m_copy.color    =grain->color;
	m_copy.selection=grain->selection;
	m_copy.slices   =grain->slices;
	m_copy.layer_distance = grain->layer_distance;
*/

	Quaternionr q(b->physicalParameters->se3.orientation);

	Vector3r min(q*(bss->m_copy.slices[0][0]));
	Vector3r max(min);
	BOOST_FOREACH(std::vector<Vector3r>& vv, bss->m_copy.slices)
	BOOST_FOREACH(Vector3r v, vv)
	{
	 	max = componentMaxVector(max,q*v);
 		min = componentMinVector(min,q*v);
	}
	
	max += b->physicalParameters->se3.position;
	min += b->physicalParameters->se3.position;

	aabb->center = (max+min)*0.5;
	aabb->halfSize = (max-min)*0.5;
	
	aabb->min = min;
	aabb->max = max;
/*
	aabb->center = se3.position;
	aabb->halfSize = 2.0*Vector3r(sphere->radius,sphere->radius,sphere->radius);
	
	aabb->min = aabb->center-aabb->halfSize;
	aabb->max = aabb->center+aabb->halfSize;	
*/
}
	
YADE_PLUGIN((Ef2_BssSnowGrain_AABB_makeAABB));

YADE_REQUIRE_FEATURE(PHYSPAR);

