/*************************************************************************
*  Copyright (C) 2006 by Sergei Dorofeenko				 *
*  sega@users.berlios.de                                                 *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include <yade/pkg-common/GLDrawFunctors.hpp>

class GLDrawViscoelasticInteraction : public GlInteractionPhysicsFunctor
{
	public :
		virtual void go(const shared_ptr<InteractionPhysics>&,const shared_ptr<Interaction>&,const shared_ptr<Body>&,const shared_ptr<Body>&,bool wireFrame);

	RENDERS(ViscoelasticInteraction);
	REGISTER_CLASS_NAME(GLDrawViscoelasticInteraction);
	REGISTER_BASE_CLASS_NAME(GlInteractionPhysicsFunctor);
};

REGISTER_SERIALIZABLE(GLDrawViscoelasticInteraction);


