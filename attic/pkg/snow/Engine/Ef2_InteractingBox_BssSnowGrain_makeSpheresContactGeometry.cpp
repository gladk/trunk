/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include"Ef2_InteractingBox_BssSnowGrain_makeSpheresContactGeometry.hpp"
#include<yade/pkg-dem/ScGeom.hpp>
#include<yade/pkg-common/Sphere.hpp>
#include<yade/pkg-common/Box.hpp>
#include<yade/pkg-snow/BssSnowGrain.hpp>

#include<yade/lib-base/yadeWm3Extra.hpp>

YADE_REQUIRE_FEATURE(geometricalmodel);

bool Ef2_InteractingBox_BssSnowGrain_makeSpheresContactGeometry::go(
		const shared_ptr<Shape>& cm1,
		const shared_ptr<Shape>& cm2,
		const Se3r& se31,
		const Se3r& se32,
		const shared_ptr<Interaction>& c)
{
	bool result = g.go(cm1,cm2,se31,se32,c) || assist;
	//std::cerr << "------------------- " << __FILE__ << "\n";

	if(result)
	{
		//Box* s1=static_cast<Box*>(cm1.get()), *s2=static_cast<BssSnowGrain*>(cm2.get());
		if(cm1->getClassName() != std::string("Box") || cm2->getClassName() != std::string("BssSnowGrain"))
		{
			std::cerr << cm1->getClassName() << " " << cm2->getClassName() << "\n";
			std::cerr << "whooooooooops =22=\n";
			return false;
		}
				Box* s1=static_cast<Box*>(cm1.get());
				BssSnowGrain *s2=static_cast<BssSnowGrain*>(cm2.get());
				if(s1==0 || s2==0)
				{
					std::cerr << cm1->getClassName() << " " << cm2->getClassName() << "\n";
					std::cerr << "whooooooooops =2= " << __FILE__ << "\n"; 
					return false;
				}
			
		shared_ptr<ScGeom> scm;
		if(c->interactionGeometry) scm=dynamic_pointer_cast<ScGeom>(c->interactionGeometry);
		else 
		{ 
			std::cerr << "whooooooooops =3= " << __FILE__ << "\n";
		}

//	std::cerr << __FILE__ << " " << scm->getClassName() << "\n";
		
		int id1 = c->getId1();
//		int id2 = c->getId2();

//		if(s1->sphere_depth.find(id2) == s1->sphere_depth.end())
//			s1->sphere_depth[id2] = scm->penetrationDepth;
		if(s2->sphere_depth.find(id1) == s2->sphere_depth.end())
			s2->sphere_depth[id1] = scm->penetrationDepth;
		
//		Real d1 = s1->sphere_depth[id2];
		Real d2 = s2->sphere_depth[id1];
//		if(d1 != d2)
//			std::cerr << "bad initial penetration?\n";
		
		scm->penetrationDepth -= d2;
	}

	return result;
}


bool Ef2_InteractingBox_BssSnowGrain_makeSpheresContactGeometry::goReverse(	const shared_ptr<Shape>& cm1,
						const shared_ptr<Shape>& cm2,
						const Se3r& se31,
						const Se3r& se32,
						const shared_ptr<Interaction>& c)
{
//	bool result = g.goReverse(cm1,cm2,se31,se32,c);
//	std::cerr << "-------------------2\n";
//	return result;
	
	bool result = go(cm2,cm1,se32,se31,c);
	if(result)
	{
		shared_ptr<ScGeom> scm;
		if(c->interactionGeometry) scm=dynamic_pointer_cast<ScGeom>(c->interactionGeometry);
		else { std::cerr << "whooooooooops =4= " << __FILE__ << "\n"; return false; }
		scm->normal *= -1.0;
		std::swap(scm->radius1,scm->radius2);
	}
	return result;
}

YADE_PLUGIN((Ef2_InteractingBox_BssSnowGrain_makeSpheresContactGeometry));

YADE_REQUIRE_FEATURE(PHYSPAR);

