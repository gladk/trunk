/**************************************************************************
*  Copyright (C) 2008 by Vincent Richefeu                                 *
*  vincent.richefeu@hmg.inpg.fr                                           *
*                                                                         *
*  This program is free software; it is licensed under the terms of the   *
*  GNU General Public License v[1] or later. See file LICENSE for details *
***************************************************************************/

#include "GLDrawTube.hpp"
#include<yade/pkg-common/BshTube.hpp>
#include<yade/lib-opengl/OpenGLWrapper.hpp>
#include<yade/lib-base/yadeWm3.hpp>

YADE_REQUIRE_FEATURE(geometricalmodel);

bool GLDrawTube::first           = true;
int  GLDrawTube::glWiredTubeList = -1;
int  GLDrawTube::glTubeList      = -1;
    
void GLDrawTube::go(const shared_ptr<GeometricalModel>& gm, const shared_ptr<PhysicalParameters>&,bool wire)
{
	
	if (first)
	{
                const Real my2PI  = 6.28318530717958647692;
		const Real myPI_6 = 5.23598775598298873077e-1;
		
		glWiredTubeList = glGenLists(1);
		glNewList(glWiredTubeList,GL_COMPILE);
			glBegin(GL_LINES);
			glDisable(GL_LIGHTING);
			for (Real a = 0.0 ; a < my2PI ; a += myPI_6)
			{
				glVertex3f(-1.0f, cos(a), sin(a));
				glVertex3f( 1.0f, cos(a), sin(a));
			}
			glEnd();
		glEndList();

		glTubeList = glGenLists(1);
		glNewList(glTubeList,GL_COMPILE);
			glDisable(GL_CULL_FACE);
			glEnable(GL_LIGHTING);
			glBegin(GL_TRIANGLE_STRIP);
			for (Real a = 0.0 ; a < my2PI ; a += myPI_6)
			{
				glVertex3f(-1.0f, cos(a), sin(a));
				glVertex3f( 1.0f, cos(a), sin(a));
			}
                        glVertex3f(-1.0f, cos(0.0), sin(0.0));
                        glVertex3f( 1.0f, cos(0.0), sin(0.0));
			glEnd();
		glEndList();

		first = false;
	}

	glMaterialv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, Vector3f(gm->diffuseColor[0],gm->diffuseColor[1],gm->diffuseColor[2]));
	glColor3v(gm->diffuseColor);

	BshTube* tube = static_cast<BshTube*>(gm.get());

	Real radius = tube->radius;
	Real half_height = tube->half_height;


	if (gm->wire || wire)
	{
		glPushMatrix();
		glScalef(half_height,radius,radius);
                glCallList(glWiredTubeList);
		glPopMatrix();
	}
	else
	{
		glPushMatrix();
		glScalef(half_height,radius,radius);
                
                glCallList(glTubeList);
		glPopMatrix();

	}
}

YADE_PLUGIN((GLDrawTube));
YADE_REQUIRE_FEATURE(OPENGL)

YADE_REQUIRE_FEATURE(PHYSPAR);

