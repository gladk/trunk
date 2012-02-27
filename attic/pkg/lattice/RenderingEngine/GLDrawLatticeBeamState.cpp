/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "GLDrawLatticeBeamState.hpp"
#include<yade/pkg-lattice/LatticeBeamParameters.hpp>
#include<yade/lib-opengl/OpenGLWrapper.hpp>

void GLDrawLatticeBeamState::go(const shared_ptr<PhysicalParameters>& pp)
	// FIXME:  instead of one argument 'bool wire_frame' it should get a class containing whole configuration of OpenGLRenderingEngine
{
	// FIXME the same in  GeometricalModelForceColorizer.cpp 
	static Real maxTensileFactor = 0.0; // FIXME - thread unsafe
	static Real maxCompressFactor = 0.0; // FIXME - thread unsafe

	static unsigned int cccc=0;
	if(((++cccc)%100000) == 0) std::cerr << maxTensileFactor << " " <<maxCompressFactor << "\n";

	LatticeBeamParameters* beam = static_cast<LatticeBeamParameters*>(pp.get());
	Real strain                     = beam->strain();
	Real factor;
	if( strain > 0 )
		factor                  = strain / beam->criticalTensileStrain; // positive
	else
		factor 			= strain / beam->criticalCompressiveStrain; // negative

	// compute optimal red/blue colors	
	maxTensileFactor		= std::max(factor,maxTensileFactor);
	maxCompressFactor 		= std::min(factor,maxCompressFactor);
	if(factor > 0 && maxTensileFactor > 0)
	{
		factor                  /= maxTensileFactor;
		glColor3v(Vector3r(0.9,0.9,1.0) - factor * Vector3r(0.9,0.9,0.0));
	}
	else if (factor < 0 && maxCompressFactor < 0)
	{
		factor                  /= maxCompressFactor; // both are negative, so result is positive
		glColor3v(Vector3r(1.0,0.9,0.9) - factor * Vector3r(0.0,0.9,0.9));
	}
	else
		glColor3(0.9,0.9,0.9);

	glTranslatev(beam->se3.position);

	glScale(beam->length,beam->length,beam->length);

	glDisable(GL_LIGHTING);
	glBegin(GL_LINES);
		glVertex3v(beam->direction*(-0.5));
		glVertex3v(beam->direction*0.5);
	glEnd();
	/*
	glBegin(GL_TRIANGLES);
		glColor3(1.0,1.0,1.0);
		glVertex3v(beam->direction*(-0.1));
		glVertex3v(beam->direction*0.1);
		glVertex3v(beam->otherDirection*0.2);
		
		glVertex3v(beam->direction*0.1);
		glVertex3v(beam->direction*(-0.1));
		glVertex3v(beam->otherDirection*0.2);
	glEnd();
	*/
}

YADE_PLUGIN((GLDrawLatticeBeamState));
YADE_REQUIRE_FEATURE(OPENGL)

YADE_REQUIRE_FEATURE(PHYSPAR);

