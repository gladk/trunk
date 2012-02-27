/*************************************************************************
*  Copyright (C) 2007 by Janek Kozicki                                   *
*  cosurgi@mail.berlios.de                                               *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "GLDrawCohesiveFrictionalContactInteraction.hpp"
#include<yade/pkg-dem/CohFrictPhys.hpp>
#include<yade/pkg-dem/ScGeom.hpp>
#include<yade/pkg-common/NormShearPhys.hpp>
#include<yade/pkg-dem/FrictPhys.hpp>

#include<yade/lib-opengl/OpenGLWrapper.hpp>
#include<yade/lib-QGLViewer/qglviewer.h>

GLDrawCohesiveFrictionalContactInteraction::GLDrawCohesiveFrictionalContactInteraction() : boxIndex(-1)
{
}

void GLDrawCohesiveFrictionalContactInteraction::drawAxis(const Vector3r pos,const Quaternionr q,Real size,const int old)
{
	if(old==1) { // stupid way to change axis' color
		const GLfloat ambientColor[4]	= {0.0,0.0,isReal?10.0:0.1,1.0};
		glLightModelfv(GL_LIGHT_MODEL_AMBIENT,ambientColor);
	} else if(old==2){
		const GLfloat ambientColor[4]	= {0.0,isReal?10.0:0.1,0.0,1.0};
		glLightModelfv(GL_LIGHT_MODEL_AMBIENT,ambientColor);
	} else if(old==4){
		const GLfloat ambientColor[4]	= {isReal?10.0:0.1,0.0,0.0,1.0};
		glLightModelfv(GL_LIGHT_MODEL_AMBIENT,ambientColor);
	} else {
		const GLfloat ambientColor[4]	= {isReal?10.0:0.1,isReal?10.0:0.1,0.0,1.0};
		glLightModelfv(GL_LIGHT_MODEL_AMBIENT,ambientColor);
	}

	Real matrix[16];

	glPushMatrix();
	glTranslatev(pos);
	quaterniontoGLMatrix(q,matrix);
	glMultMatrix(matrix);
	QGLViewer::drawAxis(size);
	glPopMatrix();

	const GLfloat ambientColor[4]	= {0.5,0.5,0.5,1.0};	
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT,ambientColor);
};
		
void GLDrawCohesiveFrictionalContactInteraction::drawArrow(const Vector3r from,const Vector3r to,const Vector3r color)
{
	glEnable(GL_LIGHTING);
	glColor3v(color);
	qglviewer::Vec a(from[0],from[1],from[2]),b(to[0],to[1],to[2]);
	QGLViewer::drawArrow(a,b);	
};

void GLDrawCohesiveFrictionalContactInteraction::drawFlatText(const Vector3r pos,const std::string txt)
{
	glPushMatrix();
	glTranslatev(pos);
	glColor3(1.0,1.0,0.0);
	glRasterPos2i(0,0);
	for(unsigned int i=0;i<txt.length();i++)
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, txt[i]);
	glPopMatrix();
};

Vector3r GLDrawCohesiveFrictionalContactInteraction::drawRotationAxis(
	const Quaternionr q,
	const std::string txt,
	const Vector3r pos,
	const Real size,
	const Vector3r color)
{
	glPushMatrix();
	Vector3r axis;Real angle;
	Quaternionr delta(q);
	delta.ToAxisAngle(axis,angle);
	if(angle==0) return Vector3r(0,0,0);
	if(angle > Mathr::PI) angle -= Mathr::TWO_PI; // angle is between 0 and 2*pi, but should be between -pi and pi 
	//Vector3r m(pos+axis*size*(angle>0?1.0:-1.0));
	Vector3r ret = axis*angle;
	Vector3r m(pos+ret*size);
	drawArrow(pos, m ,color);
	drawFlatText(m,txt);
	glPopMatrix();
	return ret;
}

void GLDrawCohesiveFrictionalContactInteraction::go(
		const shared_ptr<InteractionPhysics>& ih,
		const shared_ptr<Interaction>& i,
		const shared_ptr<Body>& b1,
		const shared_ptr<Body>& b2,
		bool wireFrame)
{
	if(!i->isReal()) return;
	isReal=i->isReal();

	CohFrictPhys*    ph = static_cast<CohFrictPhys*>(ih.get());
	ScGeom*    sc = static_cast<ScGeom*>(i->interactionGeometry.get());
	Vector3r pos1   = b1->physicalParameters->se3.position;
	Vector3r pos2   = b2->physicalParameters->se3.position;

	/// skip drawing interaction with boxes....
		if(boxIndex == -1 && b1->shape->getClassName() == "Box") boxIndex = b1->shape->getClassIndex();
		if(boxIndex == -1 && b2->shape->getClassName() == "Box") boxIndex = b2->shape->getClassIndex();
		if(b1->shape->getClassIndex() == boxIndex || b2->shape->getClassIndex() == boxIndex) return;
	//

	Vector3r middle = 0.5*(pos1+pos2);
	Vector3r dist   = 0.5*(pos2-pos1);
	Vector3r pos1d  = 0.9*pos1 + 0.1*pos2;
	Vector3r pos2d  = 0.1*pos1 + 0.9*pos2;
	Real     size   = dist.Length()*0.5;
	Vector3r color(isReal?1:0.4,isReal?1:0.4,isReal?1:0.4);

// draw 1st body
// 1st body initial orientation
	drawAxis(pos1,ph->initialOrientation1                 * ph->orientationToContact1,size,1);
// 1st body real orientation
	drawAxis(pos1,b1->physicalParameters->se3.orientation                            ,size,4);
// 1st body current orientation
	drawAxis(pos1,b1->physicalParameters->se3.orientation * ph->orientationToContact1,size);
                                                                                         
// draw 2st body                                                                         
// 2nd body initial orientation
	drawAxis(pos2,ph->initialOrientation2                 * ph->orientationToContact2,size,1);
// 2st body real orientation
	drawAxis(pos2,b2->physicalParameters->se3.orientation                            ,size,4);
// 2nd body current orientation
	drawAxis(pos2,b2->physicalParameters->se3.orientation * ph->orientationToContact2,size);

// draw middle
// middle initial orientation
//	drawAxis(middle,ph->initialContactOrientation,size,1);
// middle current orientation
//	drawAxis(middle,ph->currentContactOrientation,size);



	Quaternionr init1(ph->initialOrientation1                 * ph->orientationToContact1);
	Quaternionr init2(ph->initialOrientation2                 * ph->orientationToContact2);

	Quaternionr curr1(b1->physicalParameters->se3.orientation * ph->orientationToContact1);
	Quaternionr curr2(b2->physicalParameters->se3.orientation * ph->orientationToContact2);
	
	Vector3r M1 = drawRotationAxis(curr1 * init1.Conjugate() , "  M1", pos1, size*0.9,color);
	Vector3r M2 = drawRotationAxis(curr2 * init2.Conjugate() , "  M2", pos2, size*0.9,color);

//		drawRotationAxis(curr1 * init1.Conjugate() , "  M1", middle, size*0.9,color);
//		drawRotationAxis(curr2 * init2.Conjugate() , "  M2", middle, size*0.9,color);

if(wireFrame)
{
	drawRotationAxis(curr1 * curr2.Conjugate() , "  M", middle, size*50, Vector3r(isReal?1.0:0.4,0,0));

//	drawArrow(middle, middle+ph->shearDisplacement*size ,Vector3r(0,0,isReal?1:0.4));
//	drawFlatText(     middle+ph->shearDisplacement*size ,std::string("  shearDisplacement"));
	
//	drawArrow(middle        , middle+(M2-M1)*size ,Vector3r(2.0,0.4,0.4));
//	drawArrow(middle+M2*size, middle+(M2-M1)*size ,Vector3r(1.0,0.4,0.4));


	drawArrow(middle, middle+ph->shearForce*size/4000000.0 ,Vector3r(isReal?1:0.4,0,isReal?1:0.4));
	drawFlatText(     middle+ph->shearForce*size/4000000.0 ,std::string(" Fs"));

// draw normal
	drawArrow(middle, middle+sc->normal*size*0.9 ,Vector3r(0,isReal?1:0.4,0));
	drawFlatText(     middle+sc->normal*size*0.9 ,std::string(" N"));

// draw prevNormal
//	drawArrow(middle, middle+ph->prevNormal*size*0.9 ,Vector3r(isReal?1:0.4,0,0));
// draw initial normal
//	drawArrow(middle, middle+(ph->initialContactOrientation*Vector3r(1,0,0))*size*0.4 ,Vector3r(0,0,isReal?1:0.4));
// draw current normal
//	drawArrow(middle, middle+(ph->currentContactOrientation*Vector3r(1,0,0))*size*0.4 ,Vector3r(isReal?1:0.4,isReal?1:0.4,0));

//// draw initX1
//	drawArrow(pos1, pos1+ph->initX1*size*0.4 ,Vector3r(0,0,isReal?1:0.4));
//	drawFlatText(pos1+ph->initX1*size*0.4,std::string("  initX1"));
//// draw prevX1
//	drawArrow(pos1, pos1+ph->prevX1*size*0.4 ,Vector3r(isReal?1:0.4,isReal?1:0.4,0));
//	drawFlatText(pos1+ph->prevX1*size*0.4,std::string("  prevX1"));
//// draw initX2
//	drawArrow(pos2, pos2+ph->initX2*size*0.4 ,Vector3r(0,0,isReal?1:0.4));
//	drawFlatText(pos2+ph->initX2*size*0.4,std::string("  initX2"));
//// draw prevX2
//	drawArrow(pos2, pos2+ph->prevX2*size*0.4 ,Vector3r(isReal?1:0.4,isReal?1:0.4,0));
//	drawFlatText(pos2+ph->prevX2*size*0.4,std::string("  prevX2"));

// draw moment
//	drawRotationAxis(ph->currentContactOrientation             * ph->initialContactOrientation.Conjugate() , "  ref", middle, size*0.9,color);
}

//{// draw contact vector
//	//Quaternionr curr1(de1->se3.orientation * currentContactPhysics->orientationToContact1);
//	//Quaternionr curr2(de2->se3.orientation * currentContactPhysics->orientationToContact2);
//	Vector3r initialContactPoint1(pos1 + curr1*Vector3r(     sc->radius1,0,0) );
//	Vector3r initialContactPoint2(pos2 + curr2*Vector3r(-1.0*sc->radius2,0,0) );
//	Vector3r contactDist( initialContactPoint2 - initialContactPoint1 );
////	Vector3r shearDir(ph->shearForce); shearDir.Normalize();
////	Vector3r contactDir(contactDist); contactDir.Normalize();
////	Vector3r contactShear( shearDir * shearDir.Dot(contactDir) * contactDist.Length() );
//	drawArrow(initialContactPoint1, initialContactPoint1+contactDist ,Vector3r(0.2,0.2,0.2));
////	drawArrow(initialContactPoint1, initialContactPoint1+contactShear ,Vector3r(0.2,0.2,0.2));
//
////if((b1->getId()==22 && b2->getId()==6)){
////std::cerr << "contactShearDisplacementGL: " << contactShear << " length:" << contactShear.Length() << "\n";
////}
//
//}
// draw initial contact point 1st body
	glColor3(1.0,1.0,0.0);
	Real matrix[16];
	glPushMatrix();
	glTranslatev(pos1);
	quaterniontoGLMatrix(b1->physicalParameters->se3.orientation * ph->orientationToContact1,matrix);
	glMultMatrix(matrix);
	glTranslate(sc->radius1,(Real)0.0,(Real)0.0);
	glutSolidSphere(size*0.1,4,4);
	glPopMatrix();

// draw initial contact point 2nd body
	glPushMatrix();
	glTranslatev(pos2);
	quaterniontoGLMatrix(b2->physicalParameters->se3.orientation * ph->orientationToContact2,matrix);
	glMultMatrix(matrix);
	glTranslate(-1.0*sc->radius2,(Real)0.0,(Real)0.0);
	glutSolidSphere(size*0.1,4,4);
	glPopMatrix();

// draw connecting line.
	glBegin(GL_LINES);
	glColor3(1.4,1.4,0.4);
	glVertex3v(pos1);
	glVertex3v(pos2);
	glEnd();

// write A,B
	drawFlatText(pos1d,std::string("  A ") + boost::lexical_cast<std::string>(b1->getId()));
	drawFlatText(pos2d,std::string("  B ") + boost::lexical_cast<std::string>(b2->getId()));

}

YADE_PLUGIN((GLDrawCohesiveFrictionalContactInteraction));
YADE_REQUIRE_FEATURE(OPENGL)

YADE_REQUIRE_FEATURE(PHYSPAR);

