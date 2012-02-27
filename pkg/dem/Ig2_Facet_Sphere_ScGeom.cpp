/*************************************************************************
*  Copyright (C) 2008 by Sergei Dorofeenko				 *
*  sega@users.berlios.de                                                 *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include"Ig2_Facet_Sphere_ScGeom.hpp"
#include<yade/pkg-dem/ScGeom.hpp>
#include<yade/pkg-common/Sphere.hpp>
#include<yade/pkg-common/Facet.hpp>
#include<yade/pkg-common/Wall.hpp>
#include<yade/core/Scene.hpp>
#include<yade/lib-base/Math.hpp>

YADE_PLUGIN((Ig2_Facet_Sphere_ScGeom)(Ig2_Wall_Sphere_ScGeom));

CREATE_LOGGER(Ig2_Facet_Sphere_ScGeom);

bool Ig2_Facet_Sphere_ScGeom::go(const shared_ptr<Shape>& cm1,
							const shared_ptr<Shape>& cm2,
							const State& state1,
							const State& state2,
							const Vector3r& shift2,
							const bool& force,
							const shared_ptr<Interaction>& c)
{
	const Se3r& se31=state1.se3; const Se3r& se32=state2.se3;
	Facet*   facet = static_cast<Facet*>(cm1.get());
	/* could be written as (needs to be tested):
	 * Vector3r cl=se31.orientation.Conjugate()*(se32.position-se31.position);
	 */
	Matrix3r facetAxisT=se31.orientation.toRotationMatrix(); 
	Matrix3r facetAxis = facetAxisT.transpose();
	// local orientation
	Vector3r cl = facetAxis*(se32.position + shift2 - se31.position);  // "contact line" in facet-local coords

	//
	// BEGIN everything in facet-local coordinates
	//

	Vector3r normal = facet->nf;
	Real L = normal.dot(cl);
	if (L<0) {normal=-normal; L=-L; }

	Real sphereRadius = static_cast<Sphere*>(cm2.get())->radius;
	if (L>sphereRadius && !c->isReal() && !force)  return false; // no contact, but only if there was no previous contact; ortherwise, the constitutive law is responsible for setting Interaction::isReal=false

	Vector3r cp = cl - L*normal;
	const Vector3r* ne = facet->ne;
	
	Real penetrationDepth=0;
	
	Real bm = ne[0].dot(cp);
	int m=0;
	for (int i=1; i<3; ++i)
	{
		Real b=ne[i].dot(cp);
		if (bm<b) {bm=b; m=i;}
	}

	Real sh = sphereRadius*shrinkFactor;
	Real icr = facet->icr-sh;

	if (icr<0)
	{
		LOG_WARN("a radius of a facet's inscribed circle less than zero! So, shrinkFactor is too large and would be reduced to zero.");
		shrinkFactor=0;
		icr = facet->icr;
		sh = 0;
	}


	if (bm<icr) // contact with facet's surface
	{
		penetrationDepth = sphereRadius - L;
		normal.normalize();
	}
	else
	{
		cp = cp + ne[m]*(icr-bm);
		if (cp.dot(ne[(m-1<0)?2:m-1])>icr) // contact with vertex m
//			cp = facet->vertices[m];
			cp = facet->vu[m]*(facet->vl[m]-sh);
		else if (cp.dot(ne[m=(m+1>2)?0:m+1])>icr) // contact with vertex m+1
//			cp = facet->vertices[(m+1>2)?0:m+1];
			cp = facet->vu[m]*(facet->vl[m]-sh);
		normal = cl-cp;
		Real norm=normal.norm(); normal/=norm;
		penetrationDepth = sphereRadius - norm;
	}

	//
	// END everything in facet-local coordinates
	//

	if (penetrationDepth>0 || c->isReal())
	{
		shared_ptr<ScGeom> scm;
		bool isNew = !c->geom;
		if (c->geom)
			scm = YADE_PTR_CAST<ScGeom>(c->geom);
		else
			scm = shared_ptr<ScGeom>(new ScGeom());
	  
		normal = facetAxisT*normal; // in global orientation
		scm->contactPoint = se32.position + shift2 - (sphereRadius-0.5*penetrationDepth)*normal;
		scm->penetrationDepth = penetrationDepth;
		scm->radius1 = 2*sphereRadius;
		scm->radius2 = sphereRadius;
		if (isNew) c->geom = scm;
		scm->precompute(state1,state2,scene,c,normal,isNew,shift2,false/*avoidGranularRatcheting only for sphere-sphere*/);
		return true;
	}
	return false;
}


bool Ig2_Facet_Sphere_ScGeom::goReverse(	const shared_ptr<Shape>& cm1,
								const shared_ptr<Shape>& cm2,
								const State& state1,
								const State& state2,
								const Vector3r& shift2,
								const bool& force,
								const shared_ptr<Interaction>& c)
{
	c->swapOrder();
	//LOG_WARN("Swapped interaction order for "<<c->getId2()<<"&"<<c->getId1());
	return go(cm2,cm1,state2,state1,-shift2,force,c);
}


/********* Wall + Sphere **********/

bool Ig2_Wall_Sphere_ScGeom::go(const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c){
	Wall* wall=static_cast<Wall*>(cm1.get());
	const Real radius=static_cast<Sphere*>(cm2.get())->radius;
	const int& ax(wall->axis);
	Real dist=(state2.pos)[ax]+shift2[ax]-state1.pos[ax]; // signed "distance" between centers
	if(!c->isReal() && abs(dist)>radius && !force) { return false; }// wall and sphere too far from each other

	// contact point is sphere center projected onto the wall
	Vector3r contPt=state2.pos+shift2; contPt[ax]=state1.pos[ax];
	Vector3r normal(0.,0.,0.);
	// wall interacting from both sides: normal depends on sphere's position
	assert(wall->sense==-1 || wall->sense==0 || wall->sense==1);
	if(wall->sense==0) normal[ax]=dist>0?1.:-1.;
	else normal[ax]=wall->sense==1?1.:-1;

	bool isNew=!c->geom;
	if(isNew) c->geom=shared_ptr<ScGeom>(new ScGeom());
	const shared_ptr<ScGeom>& ws=YADE_PTR_CAST<ScGeom>(c->geom);
	ws->radius1=ws->radius2=radius; // do the same as for facet-sphere: wall's "radius" is the same as the sphere's radius
	ws->contactPoint=contPt;
	ws->penetrationDepth=-(abs(dist)-radius);
	// ws->normal is assigned by precompute
	ws->precompute(state1,state2,scene,c,normal,isNew,shift2,noRatch);
	return true;
}

