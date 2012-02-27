/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/pkg-common/InteractionGeometryFunctor.hpp>
//#include<yade/pkg-dem/Ig2_Box_Sphere_ScGeom.hpp>
#include<yade/pkg-snow/Ef2_BssSnowGrain_BssSnowGrain_makeIstSnowLayersContact.hpp>
#include<yade/pkg-snow/Ef2_InteractingBox_BssSnowGrain_makeSpheresContactGeometry.hpp>

class Ef2_InteractingBox_BssSnowGrain_makeIstSnowLayersContact : public InteractionGeometryFunctor
{
	Ef2_BssSnowGrain_BssSnowGrain_makeIstSnowLayersContact g;
	Ef2_InteractingBox_BssSnowGrain_makeSpheresContactGeometry ggg;
	public :
		virtual bool go(	const shared_ptr<Shape>& cm1,
					const shared_ptr<Shape>& cm2,
					const State& state1, const State state2, const Vector3r& shift2, const bool& force,
					const shared_ptr<Interaction>& c);

		virtual bool goReverse(	const shared_ptr<Shape>& cm1,
					const shared_ptr<Shape>& cm2,
					const State& state1, const State state2, const Vector3r& shift2, const bool& force,
					const shared_ptr<Interaction>& c);

	DECLARE_LOGGER;

	REGISTER_CLASS_NAME(Ef2_InteractingBox_BssSnowGrain_makeIstSnowLayersContact);
	REGISTER_BASE_CLASS_NAME(InteractionGeometryFunctor);

	FUNCTOR2D(Box,BssSnowGrain);

	DEFINE_FUNCTOR_ORDER_2D(Box,BssSnowGrain);
};

REGISTER_SERIALIZABLE(Ef2_InteractingBox_BssSnowGrain_makeIstSnowLayersContact);


