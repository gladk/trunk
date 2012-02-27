/*************************************************************************
*  Copyright (C) 2008 by Vincent Richefeu                                *
*  vincent.richefeu@hmg.inpg.fr                                          *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "GLDrawBssSweptSphereLineSegment.hpp"
#include<yade/pkg-common/BssSweptSphereLineSegment.hpp>
#include<yade/lib-opengl/OpenGLWrapper.hpp>

bool GLDrawBssSweptSphereLineSegment::first       = true;
int  GLDrawBssSweptSphereLineSegment::glList      = -1;
int  GLDrawBssSweptSphereLineSegment::glWiredList = -1;

GLDrawBssSweptSphereLineSegment::GLDrawBssSweptSphereLineSegment(){ /*first=true;*/ };

void GLDrawBssSweptSphereLineSegment::go(const shared_ptr<Shape>& cm, const shared_ptr<PhysicalParameters>& ,bool,const GLViewInfo&)
{
        
  if (first)
  {        
        const Real my2PI  = 6.28318530717958647692;
        const Real myPI_6 = 5.23598775598298873077e-1;
                                
        glList = glGenLists(1);
        glNewList(glList,GL_COMPILE);
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

  BssSweptSphereLineSegment* ssls = static_cast<BssSweptSphereLineSegment*> (cm.get());
  
  Real radius = ssls->radius;
  Real half_length = 0.5 * ssls->length;

  glMaterialv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, Vector3f(cm->diffuseColor[0],cm->diffuseColor[1],cm->diffuseColor[2]));
  glColor3v(cm->diffuseColor);
  
  static GLUquadric* quadric = gluNewQuadric();
  
  glPushMatrix();
  glScalef(half_length,radius,radius);
  glCallList(glList);
  glPopMatrix();
  
  // It is not very elegant but it provides assistance
  // TODO - write a gl-list that draw 2 halh spheres 
  glTranslate((double)(-half_length), 0.0, 0.0);
  gluSphere(quadric, radius, 5, 5);
  glTranslate((double)(2.0*half_length), 0.0, 0.0);
  gluSphere(quadric, radius, 5, 5);
}





YADE_PLUGIN((GLDrawBssSweptSphereLineSegment));
YADE_REQUIRE_FEATURE(OPENGL)

YADE_REQUIRE_FEATURE(PHYSPAR);

