/*************************************************************************
*  Copyright (C) 2008 by Janek Kozicki                                   *
*  cosurgi@mail.berlios.de                                               *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#pragma once

#include<yade/pkg-common/GLDrawFunctors.hpp>

class GLDrawCohesiveFrictionalContactInteraction : public GlInteractionPhysicsFunctor
{	
	int boxIndex;
	bool isReal;
	public :
		GLDrawCohesiveFrictionalContactInteraction();
		virtual void go(const shared_ptr<InteractionPhysics>&,const shared_ptr<Interaction>&,const shared_ptr<Body>&,const shared_ptr<Body>&,bool wireFrame);
		void drawAxis(const Vector3r pos,const Quaternionr q,const Real size,const int old=0);
		void drawArrow(const Vector3r from,const Vector3r to,const Vector3r color);
		void drawFlatText(const Vector3r pos,const std::string txt);
		Vector3r drawRotationAxis(const Quaternionr q,const std::string txt,const Vector3r pos,const Real size,const Vector3r color);

	RENDERS(CohFrictPhys);
	REGISTER_CLASS_NAME(GLDrawCohesiveFrictionalContactInteraction);
	REGISTER_BASE_CLASS_NAME(GlInteractionPhysicsFunctor);
};

REGISTER_SERIALIZABLE(GLDrawCohesiveFrictionalContactInteraction);

