/*************************************************************************
*  Copyright (C) 2008 by Janek Kozicki                                   *
*  cosurgi@berlios.de                                                    *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

#include"Ef1_BssSnowGrain_glDraw.hpp"
#include<yade/pkg-snow/BssSnowGrain.hpp>
#include<yade/pkg-snow/BshSnowGrain.hpp>
#include<yade/lib-opengl/OpenGLWrapper.hpp>
#include<yade/core/Omega.hpp>

YADE_REQUIRE_FEATURE(geometricalmodel);

void triangle(Vector3r a,Vector3r b, Vector3r c,Vector3r n)
{
	glNormal3v(n);
	glVertex3v(a);
	glVertex3v(b);
	glVertex3v(c);
}

void quad_bsh(Vector3r a,Vector3r b, Vector3r c, Vector3r d,Vector3r n)
{
	glNormal3v(n);
	glVertex3v(a);
	glVertex3v(b);
	glVertex3v(c);
	glVertex3v(d);
}

void Ef1_BssSnowGrain_glDraw::go(const shared_ptr<Shape>& cm, const shared_ptr<PhysicalParameters>& pp,bool wire, const GLViewInfo&)
{
//	s.go(cm,pp,wire);
//	return;

	bool surface = !wire;
	BssSnowGrain* bssgr = static_cast<BssSnowGrain*>(cm.get());
	BshSnowGrain& gr = bssgr->m_copy;

  	glMaterialv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, Vector3f(cm->diffuseColor[0],cm->diffuseColor[1],cm->diffuseColor[2]));
	glColor3v(cm->diffuseColor);

	std::vector<boost::tuple<Vector3r,Vector3r,Vector3r,Vector3r> > f(gr.get_faces_copy());
	if(surface)
	{
		glDisable(GL_CULL_FACE);
		glEnable(GL_LIGHTING);
		glShadeModel(GL_FLAT);
		glBegin(GL_TRIANGLES);
			for(size_t i = 0; i < f.size() ; ++i)
			{
				Vector3r a(get<0>(f[i]));
				Vector3r b(get<1>(f[i]));
				Vector3r c(get<2>(f[i]));
				Vector3r n(get<3>(f[i]));
				// plot the triangular face
				triangle(a,b,c,n);
				
				// plot the depth tetrahedron
				Real depth = gr.depth(i);
				Vector3r Z((a+b+c)/3.0 + n*depth);
				Vector3r N1((Z - a).Cross(a - b));
				Vector3r N2((Z - c).Cross(c - a));
				Vector3r N3((Z - b).Cross(b - c));
				
				triangle(b,a,Z,N1);
				triangle(a,c,Z,N2);
				triangle(c,b,Z,N3);
				
				// plot the parallelepiped top-face at the half depth
				Real depth2 = depth*0.5;
				Vector3r N(n*depth2);
				Vector3r A(a+N);
				Vector3r B(b+N);
				Vector3r C(c+N);
				triangle(A,B,C,n);
			}
		glEnd();
		glBegin(GL_QUADS);
			for(size_t i = 0; i < f.size() ; ++i)
			{
				Vector3r a(get<0>(f[i]));
				Vector3r b(get<1>(f[i]));
				Vector3r c(get<2>(f[i]));
				Vector3r n(get<3>(f[i]));
				// plot the parallelepiped side faces (quads)
				Real depth = gr.depth(i)*0.5;
				Vector3r N(n*depth);
				Vector3r A(a+N);
				Vector3r B(b+N);
				Vector3r C(c+N);

				Vector3r Z((a+b+c)/3.0 + n*depth);
				Vector3r N1((A - a).Cross(a - b));
				Vector3r N2((C - c).Cross(c - a));
				Vector3r N3((B - b).Cross(b - c));
				
				quad_bsh(b,a,A,B,N1);
				quad_bsh(a,c,C,A,N2);
				quad_bsh(c,b,B,C,N3);
			}
		glEnd();
		glShadeModel(GL_SMOOTH);
		glEnable(GL_CULL_FACE);
	}
	else
	{
		glDisable(GL_LIGHTING);
			glColor3v(cm->diffuseColor);
			for(size_t i=0;i < gr.slices.size();++i)
			{
				glBegin(GL_LINE_STRIP);
					for(size_t j=0 ; j < gr.slices[i].size() ; ++j)
						glVertex3v(gr.slices[i][j]);
					glVertex3v(gr.slices[i][0]);
				glEnd();
			}
		/*
		glBegin(GL_LINES);
			for(size_t i = 0; i < f.size() ; ++i)
			{
				Vector3r a(get<0>(f[i]));
				Vector3r b(get<1>(f[i]));
				Vector3r c(get<2>(f[i]));
				Vector3r n(get<3>(f[i]));
				// plot the depth tetrahedron
				Real depth = gr.depth(i);
				Vector3r Z((a+b+c)/3.0 + n*depth);
				glVertex3v(a);glVertex3v(Z);
				glVertex3v(b);glVertex3v(Z);
				glVertex3v(c);glVertex3v(Z);
				
				// plot the parallelepiped top-face at the half depth
				Real depth2 = depth*0.5;
				Vector3r N(n*depth2);
				Vector3r A(a+N);
				Vector3r B(b+N);
				Vector3r C(c+N);

				glVertex3v(A);glVertex3v(B);
				glVertex3v(B);glVertex3v(C);
				glVertex3v(C);glVertex3v(A);

				glVertex3v(a);glVertex3v(A);
				glVertex3v(b);glVertex3v(B);
				glVertex3v(c);glVertex3v(C);
			}
		glEnd();
		*/
		glEnable(GL_LIGHTING);
	}

	// draw m_quick_lookup
	size_t level = Omega::instance().isoSec;
	const std::vector<std::vector<std::vector<std::set<int> > > >& lookup(gr.m_quick_lookup);
	Vector3r min(gr.m_min),max(gr.m_max),dist(gr.m_dist);
	typedef std::set<int> t_set;

	for(size_t x=0;x<lookup.size();x++)
	for(size_t y=0;y<lookup[0].size();y++)
	for(size_t z=0;z<lookup[0][0].size();z++)
	{
		const t_set& s(lookup[x][y][z]);
		if(s.size() >= level)
		{
			glPushMatrix();
			glTranslatev(min + Vector3r(dist[0]*x , dist[1]*y, dist[2]*z ) + dist*0.5);
			glScalev(dist);
			glColor3(0.9,0.9,0.9);
			if(wire)
				glutWireCube(1.0);
			else
				glutSolidCube(1.0);
			glPopMatrix();
		}
	}

}


YADE_PLUGIN((Ef1_BssSnowGrain_glDraw));

YADE_REQUIRE_FEATURE(PHYSPAR);

