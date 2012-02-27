/*************************************************************************
*  Copyright (C) 2008 by Vincent Richefeu                                *
*  vincent.richefeu@hmg.inpg.fr                                          *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include"InteractingSphere2BssSweptSphereLineSegment4SpheresContactGeometry.hpp"
#include<yade/pkg-dem/ScGeom.hpp>
#include<yade/pkg-common/Sphere.hpp>
#include<yade/pkg-common/BssSweptSphereLineSegment.hpp>

#include<yade/lib-base/yadeWm3Extra.hpp>

// BUILDING WORK !!!!

bool InteractingSphere2BssSweptSphereLineSegment4SpheresContactGeometry::go(
		const shared_ptr<Shape>& cm1,
		const shared_ptr<Shape>& cm2,
		const Se3r& se31,
  		const Se3r& se32,
		const Vector3r& shift2,
		const bool& force,
		const shared_ptr<Interaction>& c)
{
        Sphere* s            = static_cast<Sphere*>(cm1.get());
        BssSweptSphereLineSegment* ssls = static_cast<BssSweptSphereLineSegment*>(cm2.get());
        
        Real SquaredHalfLength = 0.25 * ssls->length * ssls->length; // FIXME - ssls should store halfLength
        Vector3r cc = se31.position - se32.position;
        Vector3r proj = (cc.Dot(ssls->orientation)) * ssls->orientation; 
     
	if (proj.SquaredLength() < SquaredHalfLength) // Potential contact with the cylindrical part of ssls 
        {
                Vector3r ccn = cc - proj;           
                Real ccn_dist = ccn.Normalize();
                        
                Real overlap = ccn_dist - s->radius - ssls->radius;
                if (overlap <= 0.0)
                {
                        shared_ptr<ScGeom> scm;
                        if (c->interactionGeometry) scm = shared_ptr<ScGeom>(new ScGeom());
                        else scm = YADE_PTR_CAST<ScGeom>(c->interactionGeometry);
        
                        scm->contactPoint = se32.position + proj + (ssls->radius+0.5*overlap)*ccn;
                        scm->normal = -ccn;
                        scm->penetrationDepth = -overlap;
                        scm->radius1 = s->radius;
                        scm->radius2 = ssls->radius;
                        c->interactionGeometry = scm;
                        return true;
                }
                
                return false;
        }
        // else ... TODO

	return false;
}


bool InteractingSphere2BssSweptSphereLineSegment4SpheresContactGeometry::goReverse(	const shared_ptr<Shape>& cm1,
						const shared_ptr<Shape>& cm2,
						const Se3r& se31,
						const Se3r& se32,
						const Vector3r& shift2,
						const bool& force,
						const shared_ptr<Interaction>& c)
{
	bool isInteracting = go(cm2,cm1,se32,se31,c);
	if (isInteracting)
	{
		ScGeom* scm = static_cast<ScGeom*>(c->interactionGeometry.get());

                // Inverse the normal direction and swap the radii
		scm->normal = -scm->normal;
                Real tmpR  = scm->radius1;
		scm->radius1 = scm->radius2;
		scm->radius2 = tmpR;
	}
	return isInteracting;
}

YADE_PLUGIN((InteractingSphere2BssSweptSphereLineSegment4SpheresContactGeometry));

YADE_REQUIRE_FEATURE(PHYSPAR);

