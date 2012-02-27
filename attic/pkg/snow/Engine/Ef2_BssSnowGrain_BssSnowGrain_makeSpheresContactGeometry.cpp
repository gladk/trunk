// ÂŠ 2004 Olivier Galizzi <olivier.galizzi@imag.fr>
// ÂŠ 2004 Janek Kozicki <cosurgi@berlios.de>
// ÂŠ 2007 Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>
// ÂŠ 2008 VĂĄclav Ĺ milauer <eudoxos@arcig.cz>

#include"Ef2_BssSnowGrain_BssSnowGrain_makeSpheresContactGeometry.hpp"
#include<yade/pkg-dem/ScGeom.hpp>
#include<yade/pkg-snow/BssSnowGrain.hpp>
#include<yade/pkg-common/Sphere.hpp>

#include<yade/lib-base/yadeWm3Extra.hpp>
#include<yade/core/Omega.hpp>

YADE_REQUIRE_FEATURE(geometricalmodel);

bool Ef2_BssSnowGrain_BssSnowGrain_makeSpheresContactGeometry::go(	const shared_ptr<Shape>& cm1,
							const shared_ptr<Shape>& cm2,
							const Se3r& se31,
							const Se3r& se32,
							const Vector3r& shift2,
							const bool& force,
							const shared_ptr<Interaction>& c)
{
//	bool result = g.go(cm1,cm2,se31,se32,c);
//	std::cerr << "------------------- " << __FILE__ << "\n";
//	return result;

	if(cm1->getClassName() != std::string("BssSnowGrain") || cm2->getClassName() != std::string("BssSnowGrain"))
	{
		std::cerr << cm1->getClassName() << " " << cm2->getClassName() << "\n";
		std::cerr << "whooooooooops =99=\n";
		return false;
	}
	BssSnowGrain* s1=static_cast<BssSnowGrain*>(cm1.get());
	BssSnowGrain* s2=static_cast<BssSnowGrain*>(cm2.get());
	if(s1==0 || s2==0)
	{
		std::cerr << cm1->getClassName() << " " << cm2->getClassName() << "\n";
		std::cerr << "whooooooooops =9=\n";
		return false;
	}
	Vector3r normal=se32.position-se31.position;
	Real penetrationDepthSq=pow((s1->radius+s2->radius),2) - normal.SquaredLength();
	if (penetrationDepthSq>0 || c->isReal() || assist)
	{
		shared_ptr<ScGeom> scm;
		if(c->interactionGeometry) scm=dynamic_pointer_cast<ScGeom>(c->interactionGeometry);
		else { scm=shared_ptr<ScGeom>(new ScGeom()); c->interactionGeometry=scm; std::cerr << "new ScGeom\n";}
		if(scm==0)
			std::cerr << "missing scm\n";
	
//	std::cerr << __FILE__ << " " << scm->getClassName() << "\n";

		Real penetrationDepth=s1->radius+s2->radius-normal.Normalize(); /* Normalize() works in-place and returns length before normalization; from here, normal is unit vector */
		scm->contactPoint=se31.position+(s1->radius-0.5*penetrationDepth)*normal;//0.5*(pt1+pt2);
		scm->normal=normal;
		scm->penetrationDepth=penetrationDepth;

		int id1 = c->getId1();
		int id2 = c->getId2();

		if(s1->sphere_depth.find(id2) == s1->sphere_depth.end())
			s1->sphere_depth[id2] = penetrationDepth;
		if(s2->sphere_depth.find(id1) == s2->sphere_depth.end())
			s2->sphere_depth[id1] = penetrationDepth;

		Real d1 = s1->sphere_depth[id2];
		Real d2 = s2->sphere_depth[id1];
		if(d1 != d2)
			std::cerr << "bad initial penetration?\n";

		penetrationDepth -= d1;
		scm->penetrationDepth=penetrationDepth;

		scm->radius1=s1->radius;
		scm->radius2=s2->radius;
		return true;
	}
	return false;
}


bool Ef2_BssSnowGrain_BssSnowGrain_makeSpheresContactGeometry::goReverse(	const shared_ptr<Shape>& cm1,
								const shared_ptr<Shape>& cm2,
								const Se3r& se31,
								const Se3r& se32,
								const Vector3r& shift2,
								const bool& force,
								const shared_ptr<Interaction>& c)
{
	std::cerr << "---- goReverse ---- " << __FILE__ << "\n";
	bool result = go(cm2,cm1,se32,se31,c);
	if(result)
	{
		shared_ptr<ScGeom> scm;
		if(c->interactionGeometry) scm=dynamic_pointer_cast<ScGeom>(c->interactionGeometry);
		else { std::cerr << "whooooooooops =5= " << __FILE__ << "\n"; return false; }
		scm->normal *= -1.0;
		std::swap(scm->radius1,scm->radius2);
	}
	return result;
}

YADE_PLUGIN((Ef2_BssSnowGrain_BssSnowGrain_makeSpheresContactGeometry));

YADE_REQUIRE_FEATURE(PHYSPAR);

