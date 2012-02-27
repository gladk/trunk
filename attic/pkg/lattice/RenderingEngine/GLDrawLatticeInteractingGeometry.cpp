/*************************************************************************
*  Copyright (C) 2007 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include "GLDrawLatticeInteractingGeometry.hpp"
#include<yade/pkg-lattice/LatticeSetGeometry.hpp>
#include<yade/pkg-lattice/LatticeSetParameters.hpp>
#include<yade/pkg-lattice/LatticeBeamParameters.hpp>
#include<yade/pkg-lattice/LatticeNodeParameters.hpp>
#include<yade/pkg-lattice/LatticeInteractingGeometry.hpp>
#include<yade/lib-opengl/OpenGLWrapper.hpp>
#include<yade/core/Scene.hpp>
#include<yade/core/Omega.hpp>
#include<yade/lib-base/yadeWm3Extra.hpp>

YADE_REQUIRE_FEATURE(geometricalmodel);

GLDrawLatticeInteractingGeometry::GLDrawLatticeInteractingGeometry() : maxLen(0)
{
};

bool GLDrawLatticeInteractingGeometry::drawn(unsigned int A,unsigned int B,unsigned int C)
{
/*	
	unsigned int a = std::max(A,std::max(B,C));
	unsigned int c = std::min(A,std::min(B,C));
	unsigned int b;
	     if( A < a && A > c ) b=A;
	else if( B < a && B > c ) b=B;
	else if( C < a && C > c ) b=C;
	else exit(1);

//	std::cerr << a <<" " << b << " " << c << "\n";

	if(done.size()<=a) done.resize(a+1);
	if(done[a].size()<=b) done[a].resize(b+1);
	if(done[a][b].size()<=c) done[a][b].resize(c+1,0);

	if(done[a][b][c] == 1)
		return true;
	else
		done[a][b][c]=1;
	return false;
*/
	std::set<unsigned int> nodes;
	nodes.insert(A);
	nodes.insert(B);
	nodes.insert(C);
	return !done.insert(nodes).second;

};

void GLDrawLatticeInteractingGeometry::zeroDrawn()
{
	done.clear();
/*	
	std::vector<std::vector<std::vector<unsigned char> > >::iterator a=done.begin();
	std::vector<std::vector<std::vector<unsigned char> > >::iterator A=done.end();
	for( ; a!=A ; ++a)
	{
		std::vector<std::vector<unsigned char> >::iterator b=a->begin();
		std::vector<std::vector<unsigned char> >::iterator B=a->end();
		for( ; b!=B ; ++b)
		{
			std::vector<unsigned char>::iterator c=b->begin();
			std::vector<unsigned char>::iterator C=b->end();
			for( ; c!=C ; ++c)
				*c=0;
		}
	}
*/
};

void GLDrawLatticeInteractingGeometry::damagedNeighbor(unsigned int a)
{
	if(damaged.size()<=a) damaged.resize(a+1,0);
	++damaged[a];

	//(*(bodies))[A]->geometricalModel->diffuseColor = Vector3r(0.0,0.0,((float)section)/2.0);
};

void GLDrawLatticeInteractingGeometry::go(const shared_ptr<Shape>& gm, const shared_ptr<PhysicalParameters>&,bool wire, const GLViewInfo&)
{
//	return;

	static bool first=true;
	if(first)
	{
		first=false;
		Omega::instance().isoSec=1;
	}
	int limit = Omega::instance().isoSec;

	InteractionContainer* interactions = Omega::instance().getScene()->interactions.get();
	BodyContainer* bodies = Omega::instance().getScene()->bodies.get();

	InteractionContainer::iterator angles     = interactions->begin();
	InteractionContainer::iterator angles_end = interactions->end();

	if(wire)
	{
		damaged.clear();
		for(  ; angles != angles_end; ++angles )
		{
			if(        (!(bodies->exists( (*angles)->getId1())) )  
				&&    bodies->exists( (*angles)->getId2() ) )
			{
				LatticeBeamParameters* b2 = static_cast<LatticeBeamParameters*>(((*(bodies))[(*angles)->getId2()])->physicalParameters.get());
				unsigned int A = b2->id1,
					     B = b2->id2;
				Vector3r c1 = (*(bodies))[A]->geometricalModel->diffuseColor;
				Vector3r c2 = (*(bodies))[B]->geometricalModel->diffuseColor;
				// beams that have nodes with red=0 and green=0 in geometrical model are used by MovingSupport, so I want to exclude them...
				if( !((c1[0]==0 && c1[1]==0) || (c2[0]==0 && c2[1]==0)) )
					damagedNeighbor((*angles)->getId2());
			}
			else
			if(        (!(bodies->exists( (*angles)->getId2())) )  
				&&    bodies->exists( (*angles)->getId1() ) )
			{
				LatticeBeamParameters* b1 = static_cast<LatticeBeamParameters*>(((*(bodies))[(*angles)->getId1()])->physicalParameters.get());
				unsigned int A = b1->id1,
					     B = b1->id2;
				Vector3r c1 = (*(bodies))[A]->geometricalModel->diffuseColor;
				Vector3r c2 = (*(bodies))[B]->geometricalModel->diffuseColor;
				if( !((c1[0]==0 && c1[1]==0) || (c2[0]==0 && c2[1]==0)) )
					damagedNeighbor((*angles)->getId1());
			}
		}
		angles     = interactions->begin();
		angles_end = interactions->end();
	}

  	glMaterialv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, Vector3f(gm->diffuseColor[0],gm->diffuseColor[1],gm->diffuseColor[2]));
	glColor3(0.5,0.5,0.9);
	glDisable(GL_CULL_FACE);
	glEnable(GL_LIGHTING);
	glEnable(GL_NORMALIZE); // GLDrawTetrahedron
	glBegin(GL_TRIANGLES);

		zeroDrawn();
		for(  ; angles != angles_end; ++angles )
		{
			if(        bodies->exists( (*angles)->getId1() )  
				&& bodies->exists( (*angles)->getId2() ) )
			{
				bool red=false;

				LatticeBeamParameters* b1 = static_cast<LatticeBeamParameters*>(((*(bodies))[(*angles)->getId1()])->physicalParameters.get());
				LatticeBeamParameters* b2 = static_cast<LatticeBeamParameters*>(((*(bodies))[(*angles)->getId2()])->physicalParameters.get());
				maxLen = std::max(b1->length,maxLen);
				maxLen = std::max(b2->length,maxLen);
			
				Vector3r p1,p2,p3;
				unsigned int A=b1->id1,
					     B=b1->id2,
					     C=b2->id1,
					     D=b2->id2;
				if( !(A!=B && A!=C && B!=C) ){
					if( A==B )	A=D;
					else if(A==C)	A=D;
					else if(B==C)	B=D;
				}

				if(wire && (damaged.size()>static_cast<unsigned int>(std::max((*angles)->getId1(),(*angles)->getId2()))) && (damaged[(*angles)->getId1()]>=limit && damaged[(*angles)->getId2()]>=limit))
				{
					glColor3(0.9,0.2,0.2);
					red=true;
				}


				if((wire && red) || !wire)
				if(drawn(A,B,C))
				{
					p1 = static_cast<LatticeNodeParameters*>(((*(bodies))[A])->physicalParameters.get())->se3.position;
					p2 = static_cast<LatticeNodeParameters*>(((*(bodies))[B])->physicalParameters.get())->se3.position;
					p3 = static_cast<LatticeNodeParameters*>(((*(bodies))[C])->physicalParameters.get())->se3.position;
	
					Vector3r normal = (p2 - p1).UnitCross(p3 - p1);
					glNormal3v(normal);
					glVertex3v(p1+normal*maxLen*0.01);
					glVertex3v(p2+normal*maxLen*0.01);
					glVertex3v(p3+normal*maxLen*0.01);
	
					normal*=-1.0;
					glNormal3v(normal);
					glVertex3v(p1+normal*maxLen*0.01);
					glVertex3v(p2+normal*maxLen*0.01);
					glVertex3v(p3+normal*maxLen*0.01);
				}
			}
			else
			{ // deleted nearby....
			}
		}

	glEnd();

}

YADE_PLUGIN((GLDrawLatticeInteractingGeometry));
YADE_REQUIRE_FEATURE(OPENGL)

YADE_REQUIRE_FEATURE(PHYSPAR);

