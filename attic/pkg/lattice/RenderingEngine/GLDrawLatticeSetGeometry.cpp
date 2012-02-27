/*************************************************************************
*  Copyright (C) 2004 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "GLDrawLatticeSetGeometry.hpp"
#include<yade/pkg-lattice/LatticeSetGeometry.hpp>
#include<yade/pkg-lattice/LatticeSetParameters.hpp>
#include<yade/pkg-lattice/LatticeBeamParameters.hpp>
#include<yade/pkg-lattice/LatticeNodeParameters.hpp>
#include<yade/lib-opengl/OpenGLWrapper.hpp>
#include<yade/core/Scene.hpp>
#include<yade/core/Omega.hpp>
#include<yade/lib-base/yadeWm3Extra.hpp>

YADE_REQUIRE_FEATURE(geometricalmodel);

void GLDrawLatticeSetGeometry::calcMinMax()
{
/*	min=Vector3r(0,0,0);
	max=Vector3r(0.18,0.18,0.04);
	sizeX=10;
	sizeY=10;
	sizeZ=3;
*/
	min = Vector3r(1,1,1)*10000.0;
	max = min*(-1.0);
	
	Scene * lattice = static_cast<Scene*>(Omega::instance().getScene().get());
	int nodeGroupMask  = static_cast<LatticeSetParameters*>(lattice->physicalParameters.get())->nodeGroupMask;
	BodyContainer* bodies = lattice->bodies.get();
	BodyContainer::iterator bi    = bodies->begin();
	BodyContainer::iterator biEnd = bodies->end();
	for(  ; bi!=biEnd ; ++bi )  // loop over all nodes
	{
		Body* body = (*bi).get();
		if( ! ( body->getGroupMask() & nodeGroupMask ) )
			continue; // skip non-beams
		// next node
		LatticeNodeParameters* node = static_cast<LatticeNodeParameters*>(body->physicalParameters.get() );
		
	 	max = componentMaxVector(max,node->se3.position);
 		min = componentMinVector(min,node->se3.position);
	}
	float mm = std::min(max[0]-min[0],std::min(max[1]-min[1],max[2]-min[2]));
	sizeX = (int)((float)Omega::instance().isoSec*(float)(max[0]-min[0])/mm);
	sizeY = (int)((float)Omega::instance().isoSec*(float)(max[1]-min[1])/mm);
	sizeZ = (int)((float)Omega::instance().isoSec*(float)(max[2]-min[2])/mm);
	
	float dx = (max[0]-min[0])/((float)(sizeX));
	float dy = (max[1]-min[1])/((float)(sizeY));
	float dz = (max[2]-min[2])/((float)(sizeZ));

	isoStep=Vector3r(dx,dy,dz);
	max += Vector3r(dx,dy,dz)*3.0;
	min -= Vector3r(dx,dy,dz)*1.0;

	sizeX += 4;
	sizeY += 4;
	sizeZ += 4;

//	max[0] += (max[0]-min[0])/((float)(sizeX));
//	max[1] += (max[1]-min[1])/((float)(sizeX));
//	max[2] += (max[2]-min[2])/((float)(sizeX));

}

void GLDrawLatticeSetGeometry::drawPoint(Vector3r a,float color,float thickness)
{
	float I = ((a[0]-(min[0]+isoStep[0]))/((max[0]-isoStep[0]*3.0)-(min[0]+isoStep[0])))*((float)(sizeX-4))+1.5;
	float J = ((a[1]-(min[1]+isoStep[1]))/((max[1]-isoStep[1]*3.0)-(min[1]+isoStep[1])))*((float)(sizeY-4))+1.5;
	float K = ((a[2]-(min[2]+isoStep[2]))/((max[2]-isoStep[2]*3.0)-(min[2]+isoStep[2])))*((float)(sizeZ-4))+1.5;
	float t2=std::pow(thickness,2);
	int mi=std::max(0,(int)(std::floor(I)-thickness-1));int Mi=std::min(sizeX,(int)(std::ceil(I)+thickness+1));
	int mj=std::max(0,(int)(std::floor(J)-thickness-1));int Mj=std::min(sizeY,(int)(std::ceil(J)+thickness+1));
	int mk=std::max(0,(int)(std::floor(K)-thickness-1));int Mk=std::min(sizeZ,(int)(std::ceil(K)+thickness+1));
	for(int i=mi; i<Mi; i++)
		for(int j=mj; j<Mj; j++)
			for(int k=mk; k<Mk; k++)
			{
				if(   std::pow(i-I,2)+std::pow(j-J,2)+std::pow(k-K,2) < t2 )
				{
					//scalarField[i][j][k] = color;
					scalarField[i][j][k] += color; // dodawaæ a na koniec u¶redniæ
					++weights[i][j][k];
				}
			}
}

void GLDrawLatticeSetGeometry::drawLine(Vector3r a,Vector3r b,float color,float thickness)
{
//	Vector3r d = (b-a)*0.125;
//	for( float i=0 ; i<=8  ; a+=d, ++i )
//	Vector3r d = (b-a)*0.2;
//	for( float i=0 ; i<=5  ; a+=d, ++i )
	Vector3r d = (b-a)*0.333333;
	for( float i=0 ; i< 3  ; a+=d, ++i )
		drawPoint(a,color,thickness);
	
}

void GLDrawLatticeSetGeometry::generateScalarField()
{
	for(int i=0;i<sizeX;i++)
		for(int j=0;j<sizeY;j++)
			for(int k=0;k<sizeZ;k++)
			{
				scalarField[i][j][k] = 0;//(float)(i*j*k)*0.001;
				weights[i][j][k]=0;
			}


	Scene * lattice = static_cast<Scene*>(Omega::instance().getScene().get());

//	int nodeGroupMask  = static_cast<LatticeSetParameters*>(lattice->physicalParameters.get())->nodeGroupMask;
	int beamGroupMask  = static_cast<LatticeSetParameters*>(lattice->physicalParameters.get())->beamGroupMask;

	BodyContainer* bodies = lattice->bodies.get();

	BodyContainer::iterator bi    = bodies->begin();
	BodyContainer::iterator biEnd = bodies->end();

	for(  ; bi!=biEnd ; ++bi )  // loop over all beams
	{
		Body* body = (*bi).get();
		if( ! ( body->getGroupMask() & beamGroupMask ) )
			continue; // skip non-beams
		// next beam
		LatticeBeamParameters* beam = static_cast<LatticeBeamParameters*>(body->physicalParameters.get() );
		LatticeNodeParameters* node1 = static_cast<LatticeNodeParameters*>(((*(bodies))[beam->id1])->physicalParameters.get());
		LatticeNodeParameters* node2 = static_cast<LatticeNodeParameters*>(((*(bodies))[beam->id2])->physicalParameters.get());

		drawLine(node1->se3.position , node2->se3.position , beam->longitudalStiffness , Omega::instance().isoThick);
	}
	
	for(int i=0;i<sizeX;i++)
		for(int j=0;j<sizeY;j++)
			for(int k=0;k<sizeZ;k++)
				scalarField[i][j][k] /= weights[i][j][k];
}

void GLDrawLatticeSetGeometry::go(const shared_ptr<GeometricalModel>& gm, const shared_ptr<PhysicalParameters>&,bool wire)
{
	return;
	static bool initialized=false;
	if(!initialized)
	{
		Omega::instance().isoValue=1.55;
		Omega::instance().isoSec=30;
		Omega::instance().isoThick=1.0;
		oldIsoValue=0;
		oldIsoSec=0;
		oldIsoThick=0;
		initialized=true;
	}
	if(Omega::instance().isoValue!=oldIsoValue || Omega::instance().isoSec!=oldIsoSec || Omega::instance().isoThick!=oldIsoThick )
	{
		oldIsoValue=Omega::instance().isoValue;
		oldIsoSec  =Omega::instance().isoSec;
		oldIsoThick=Omega::instance().isoThick;
		calcMinMax();
		std::cerr<< "iso: " << Omega::instance().isoValue << " sect: " << Omega::instance().isoSec << " thick: " << Omega::instance().isoThick
			<< "\nmin/max:" << min << " / " << max 
			<< "\nsx/sy/sz: " << sizeX << " / " << sizeY << " / " << sizeZ << "\n";
		mc.init(sizeX,sizeY,sizeZ,min,max);
		mc.resizeScalarField(scalarField,sizeX,sizeY,sizeZ);	
		mc.resizeScalarField(weights,sizeX,sizeY,sizeZ);	
		generateScalarField();
		mc.computeTriangulation(scalarField,Omega::instance().isoValue);
	}

	// FIXME : check that : one of those 2 lines are useless
  	glMaterialv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, Vector3f(gm->diffuseColor[0],gm->diffuseColor[1],gm->diffuseColor[2]));
	glColor3v(gm->diffuseColor);
	
	//Real len = (static_cast<LatticeSetGeometry*>(gm.get()))->length;

	// FIXME - there must be a way to tell this from outside
//	glScalef(len,0.010,0.010); // it's a box, not a line. looks better :)
	//glScalef(len,0.001,0.001); // it's a box, not a line. looks better :)
/*
	if (gm->wire || wire)
	{
		glBegin(GL_LINES);
		glDisable(GL_LIGHTING);

		glVertex3(-0.5,0.0,0.0);
		glVertex3( 0.5,0.0,0.0);

		glEnd();
	}
	else
	{
		glEnable(GL_LIGHTING);
		glutSolidCube(1);
	}
	*/


	{
		const vector<Vector3r>& triangles 	= mc.getTriangles();
		int nbTriangles				= mc.getNbTriangles();
		const vector<Vector3r>& normals 	= mc.getNormals();	

		glDisable(GL_CULL_FACE);
		glEnable(GL_LIGHTING); // 2D
		//glFrontFace(GL_CCW);//: GL_CW);
		glEnable(GL_NORMALIZE);
		glBegin(GL_TRIANGLES);
			glColor3f(0.4,0.4,1.0);
			for(int i=0;i<3*nbTriangles;++i)
			{
				glNormal3v(normals[i]);
				glVertex3v(triangles[i]);
				glNormal3v(normals[++i]);
				glVertex3v(triangles[i]);
				glNormal3v(normals[++i]);
				glVertex3v(triangles[i]);
			}
		glEnd();
	
		//Vector3r size = max-min;
		//glPushMatrix();
		//glDisable(GL_LIGHTING);
		//glColor3f(1.0,1.0,1.0);
		//glScalef(size[0],size[1],size[2]);
		//glutWireCube(1);
		//glEnable(GL_LIGHTING); // 2D
		//glPopMatrix();
	}

//	std::cerr << "GLDrawLatticeSetGeometry\n";
}

YADE_PLUGIN((GLDrawLatticeSetGeometry));
YADE_REQUIRE_FEATURE(OPENGL)

YADE_REQUIRE_FEATURE(PHYSPAR);

