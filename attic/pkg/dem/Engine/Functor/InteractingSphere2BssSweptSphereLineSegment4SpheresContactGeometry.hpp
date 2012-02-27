/*************************************************************************
*  Copyright (C) 2008 by Vincent Richefeu                                *
*  vincent.richefeu@hmg.inpg.fr                                          *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/pkg-common/InteractionGeometryFunctor.hpp>

class InteractingSphere2BssSweptSphereLineSegment4SpheresContactGeometry : public InteractionGeometryFunctor
{
	public :
		virtual bool go(	const shared_ptr<Shape>& cm1,
					const shared_ptr<Shape>& cm2,
					const State& state1,
					const State& state2,
					const Vector3r& shift2,
					const bool& force,
					const shared_ptr<Interaction>& c);

		virtual bool goReverse(	const shared_ptr<Shape>& cm1,
					const shared_ptr<Shape>& cm2,
					const State& state1,
					const State& state2,
					const Vector3r& shift2,
					const bool& force,
					const shared_ptr<Interaction>& c);

        REGISTER_CLASS_NAME(InteractingSphere2BssSweptSphereLineSegment4SpheresContactGeometry);
	REGISTER_BASE_CLASS_NAME(InteractionGeometryFunctor);

        FUNCTOR2D(Sphere,BssSweptSphereLineSegment);

        DEFINE_FUNCTOR_ORDER_2D(Sphere,BssSweptSphereLineSegment);
};

REGISTER_SERIALIZABLE(InteractingSphere2BssSweptSphereLineSegment4SpheresContactGeometry);


