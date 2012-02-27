/*************************************************************************
*  Copyright (C) 2006 by Sergei Dorofeenko				 *
*  sega@users.berlios.de                                                 *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "GLDrawViscoelasticInteraction.hpp"
#include <yade/pkg-dem/ViscoelasticInteraction.hpp>
#include<yade/lib-opengl/OpenGLWrapper.hpp>


void GLDrawViscoelasticInteraction::go(
		const shared_ptr<InteractionPhysics>& ip,
		const shared_ptr<Interaction>&,
		const shared_ptr<Body>& b1,
		const shared_ptr<Body>& b2,
		bool wireFrame)
{
	ViscoelasticInteraction* p = static_cast<ViscoelasticInteraction*>(ip.get());
	Vector3r mid = 0.5*(b1->physicalParameters->se3.position+b2->physicalParameters->se3.position);
	glTranslatev(mid);
	glColor3(1.0,0.0,0.0);
	
	glPushMatrix();
	glRasterPos2i(0,0);
	std::string str;
	str = std::string("  kn=") + boost::lexical_cast<std::string>((float)(p->kn))
	    + std::string("  cn=") + boost::lexical_cast<std::string>((float)(p->cn))
	    + std::string("  ks=") + boost::lexical_cast<std::string>((float)(p->ks))
	    + std::string("  cs=") + boost::lexical_cast<std::string>((float)(p->cs));
	for(unsigned int i=0;i<str.length();i++)
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, str[i]);
	glPopMatrix();
}

YADE_PLUGIN((GLDrawViscoelasticInteraction));
YADE_REQUIRE_FEATURE(OPENGL)

YADE_REQUIRE_FEATURE(PHYSPAR);

