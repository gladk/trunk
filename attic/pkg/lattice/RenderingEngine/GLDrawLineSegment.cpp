/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "GLDrawLineSegment.hpp"
#include<yade/pkg-lattice/LineSegment.hpp>
#include<yade/lib-opengl/OpenGLWrapper.hpp>
YADE_REQUIRE_FEATURE(geometricalmodel);

GLDrawLineSegment::GLDrawLineSegment() : maxLen(0)
{
};

void GLDrawLineSegment::go(const shared_ptr<GeometricalModel>& gm, const shared_ptr<PhysicalParameters>&,bool wire)
{
	// FIXME : check that : one of those 2 lines are useless
  	glMaterialv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, Vector3f(gm->diffuseColor[0],gm->diffuseColor[1],gm->diffuseColor[2]));
	glColor3v(gm->diffuseColor);


	LineSegment* lll = dynamic_cast<LineSegment*>(gm.get());
	if(lll)
	{
	
	//Real len = (static_cast<LineSegment*>(gm.get()))->length;
	Real len = lll->length;
	//if(len >1.0)
	//	return;

	maxLen = std::max(len,maxLen);

	// FIXME - there must be a way to tell this from outside
//	glScalef(len,0.030,0.030); // it's a box, not a line. looks better :)
//	glScalef(len,0.010,0.010); // it's a box, not a line. looks better :)
//	glScalef(len,0.001,0.001); // it's a box, not a line. looks better :)
//	glScalef(len,0.00007,0.00007); // it's a box, not a line. looks better :)
	glScalef(len,maxLen*0.07,maxLen*0.07); // it's a box, not a line. looks better :)

	if (gm->wire || wire)
	{
		maxLen = 0;

		glBegin(GL_LINES);
		glDisable(GL_LIGHTING);

		glVertex3(-0.5,0.0,0.0);
		glVertex3( 0.5,0.0,0.0);

//		glVertex3( 0.0,0.0,0.0);
//		glVertex3( 0.0, 0.4*len/0.042 ,0.0);

//		glVertex3( 0.0,0.0,0.0);
//		glVertex3( 0.0, 0.0 , 0.25*len/0.042 );

		glEnd();
	}
	else
	{
		glEnable(GL_LIGHTING);
		glutSolidCube(1);
	}
	}
}

YADE_PLUGIN((GLDrawLineSegment));
YADE_REQUIRE_FEATURE(OPENGL)

YADE_REQUIRE_FEATURE(PHYSPAR);

